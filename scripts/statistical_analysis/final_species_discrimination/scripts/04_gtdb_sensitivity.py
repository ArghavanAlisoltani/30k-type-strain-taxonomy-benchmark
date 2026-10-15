#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, sys
from sklearn.metrics import roc_auc_score, average_precision_score
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'scripts'))
from common import load_config, score_for_auc, threshold_pred, confusion_metrics, genus_confusion_table, bootstrap_confusion_by_group, bootstrap_macro_mean
CFG=load_config(ROOT); IN=ROOT/'inputs/derived/binary_analysis_core.tsv'; OUT=ROOT/'outputs/tables'
df=pd.read_csv(IN,sep='\t',low_memory=False)
subsets={
 'GTDB_fixed_pair_relabeling': df[df.y_gtdb.notna()].copy(),
 'NCBI_GTDB_pair_relation_consensus': df[df.y_gtdb.notna() & df.pair_relation_consensus.eq(True)].copy(),
 'NCBI_GTDB_normalized_name_consensus': df[df.y_gtdb.notna() & df.pair_name_consensus_normalized.eq(True)].copy(),
 'NCBI_GTDB_strict_name_consensus': df[df.y_gtdb.notna() & df.pair_name_consensus_strict.eq(True)].copy(),
}
rows=[]
for si,(subset,d0) in enumerate(subsets.items()):
    for mi,(method,spec) in enumerate(CFG['metrics'].items()):
        col=spec['column']; d=d0[d0[col].notna() & d0.gtdb_genus_eval.notna()].copy(); d['y']=d.y_gtdb.astype(int)
        both=d.groupby('gtdb_genus_eval').y.nunique(); d=d[d.gtdb_genus_eval.isin(both[both==2].index)].copy()
        if d.y.nunique()<2: continue
        s=score_for_auc(d[col],spec['direction']); row={'subset':subset,'method':method,'n_pairs':len(d),'n_positive':int(d.y.sum()),'n_negative':int((1-d.y).sum()),'n_genera':d.gtdb_genus_eval.nunique(),'roc_auc_micro':roc_auc_score(d.y,s),'average_precision':average_precision_score(d.y,s)}
        # Equal-lineage macro AUC only for the main fixed-pair GTDB relabeling.
        if subset=='GTDB_fixed_pair_relabeling':
            aucs=[]
            for g,z in d.groupby('gtdb_genus_eval',sort=False):
                if int(z.y.sum())<CFG['min_lineage_positive_pairs'] or int((1-z.y).sum())<CFG['min_lineage_negative_pairs']: continue
                aucs.append(roc_auc_score(z.y,score_for_auc(z[col],spec['direction'])))
            if aucs:
                m,lo,hi,n=bootstrap_macro_mean(aucs,B=CFG['macro_auc_bootstrap_replicates'],seed=CFG['random_seed']+mi); row.update({'macro_genus_auc_mean':m,'macro_genus_auc_ci_low':lo,'macro_genus_auc_ci_high':hi,'macro_genus_n':n})
        t=spec.get('literature_threshold')
        if t is not None:
            pred=threshold_pred(d[col],t,spec['literature_rule']); d['pred']=pred; row.update({f'literature_{k}':v for k,v in confusion_metrics(d.y,pred).items()})
            # Cluster CIs only on the main fixed-pair GTDB sensitivity universe.
            if subset=='GTDB_fixed_pair_relabeling':
                gt=genus_confusion_table(d,'pred','gtdb_genus_eval','y'); ci=bootstrap_confusion_by_group(gt,B=CFG['binary_cluster_bootstrap_replicates'],seed=CFG['random_seed']+mi); row.update({f'literature_{k}':v for k,v in ci.items()})
        rows.append(row)
out=pd.DataFrame(rows); out.to_csv(OUT/'gtdb_sensitivity_summary.tsv',sep='\t',index=False)
eligible=df[df.y_gtdb.notna()].copy(); trans=pd.crosstab(eligible.y_ncbi,eligible.y_gtdb).rename_axis(index='NCBI_label',columns='GTDB_label').stack().reset_index(name='pair_count'); trans.to_csv(OUT/'NCBI_GTDB_pair_label_transition.tsv',sep='\t',index=False)
audit=pd.DataFrame([{'fixed_ncbi_pair_universe':len(df),'gtdb_binary_eligible':int(df.y_gtdb.notna().sum()),'gtdb_binary_ineligible_or_different_genus':int(df.y_gtdb.isna().sum()),'pair_relation_consensus':int((df.y_gtdb.notna()&df.pair_relation_consensus.eq(True)).sum()),'pair_relation_discordant':int((df.y_gtdb.notna()&~df.pair_relation_consensus.eq(True)).sum()),'normalized_name_consensus_pairs':int((df.y_gtdb.notna()&df.pair_name_consensus_normalized.eq(True)).sum()),'strict_name_consensus_pairs':int((df.y_gtdb.notna()&df.pair_name_consensus_strict.eq(True)).sum())}]); audit.to_csv(OUT/'gtdb_fixed_pair_sensitivity_audit.tsv',sep='\t',index=False)
print(audit.to_string(index=False)); print('\n',out[['subset','method','n_pairs','n_positive','n_negative','n_genera','roc_auc_micro','average_precision']].to_string(index=False))
