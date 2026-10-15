#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, sys
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, average_precision_score, auc
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'scripts'))
from common import load_config, score_for_auc, threshold_pred, confusion_metrics, youden_threshold, roc_pr, genus_confusion_table, bootstrap_confusion_by_group, bootstrap_macro_mean
CFG=load_config(ROOT); IN=ROOT/'inputs/derived/binary_analysis_core.tsv'; OUT=ROOT/'outputs/tables'; OUT.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(IN,sep='\t',low_memory=False).rename(columns={'y_ncbi':'y'})
summary=[]; lineage=[]
curve_dir=OUT/'curve_points'; curve_dir.mkdir(exist_ok=True)
for mi,(method,spec) in enumerate(CFG['metrics'].items()):
    col=spec['column']; d=df[df[col].notna() & df.genus_eval.notna()].copy()
    both=d.groupby('genus_eval').y.nunique(); keep=both[both==2].index; d=d[d.genus_eval.isin(keep)].copy()
    s=score_for_auc(d[col],spec['direction']); roc=float(roc_auc_score(d.y,s)); ap=float(average_precision_score(d.y,s)); p,r,_=precision_recall_curve(d.y,s); pr=float(auc(r,p)); thr_desc,j=youden_threshold(d.y,d[col],spec['direction'])
    fpr,tpr,rt=roc_curve(d.y,s); pd.DataFrame({'fpr':fpr,'tpr':tpr,'score_threshold_transformed':rt}).to_csv(curve_dir/f'{method}_roc.tsv',sep='\t',index=False)
    pd.DataFrame({'recall':r,'precision':p}).to_csv(curve_dir/f'{method}_pr.tsv',sep='\t',index=False)
    row={'method':method,'evaluation_universe':'method-specific numeric scores; NCBI genera with both classes','n_pairs':len(d),'n_positive':int(d.y.sum()),'n_negative':int((1-d.y).sum()),'n_genera':d.genus_eval.nunique(),'positive_prevalence':float(d.y.mean()),'roc_auc_micro':roc,'average_precision':ap,'pr_auc_trapezoid':pr,'full_data_youden_cutoff_descriptive_only':thr_desc,'full_data_youden_J':j}
    # genus macro AUC, deliberately equal-weighting lineages
    for g,z in d.groupby('genus_eval',sort=False):
        npos=int(z.y.sum()); nneg=int((1-z.y).sum())
        if npos<CFG['min_lineage_positive_pairs'] or nneg<CFG['min_lineage_negative_pairs']: continue
        if z.y.nunique()<2: continue
        aucg=float(roc_auc_score(z.y,score_for_auc(z[col],spec['direction']))); apg=float(average_precision_score(z.y,score_for_auc(z[col],spec['direction'])))
        lineage.append({'method':method,'genus':g,'n_positive':npos,'n_negative':nneg,'roc_auc':aucg,'average_precision':apg})
    lg=pd.DataFrame([x for x in lineage if x['method']==method])
    if len(lg):
        m,lo,hi,n=bootstrap_macro_mean(lg.roc_auc,B=CFG['macro_auc_bootstrap_replicates'],seed=CFG['random_seed']+mi); row.update({'macro_genus_auc_mean':m,'macro_genus_auc_ci_low':lo,'macro_genus_auc_ci_high':hi,'macro_genus_n':n})
    t=spec.get('literature_threshold')
    if t is not None:
        pred=threshold_pred(d[col],t,spec['literature_rule']); d['pred']=pred
        row.update({'literature_threshold':t,'literature_rule':spec['literature_rule'],**{f'literature_{k}':v for k,v in confusion_metrics(d.y,pred).items()}})
        gt=genus_confusion_table(d,'pred','genus_eval','y'); ci=bootstrap_confusion_by_group(gt,B=CFG['binary_cluster_bootstrap_replicates'],seed=CFG['random_seed']+mi)
        row.update({f'literature_{k}':v for k,v in ci.items()})
    summary.append(row)

pd.DataFrame(summary).to_csv(OUT/'primary_roc_pr_threshold_summary_ncbi.tsv',sep='\t',index=False)
lineage_df=pd.DataFrame(lineage); lineage_df.to_csv(OUT/'per_genus_auc_primary.tsv',sep='\t',index=False)
# Common complete-case comparison removes method-specific missingness as an explanation for differences.
metric_cols=[s['column'] for s in CFG['metrics'].values()]
c=df.dropna(subset=metric_cols+['genus_eval']).copy(); both=c.groupby('genus_eval').y.nunique(); c=c[c.genus_eval.isin(both[both==2].index)].copy()
cc=[]
for method,spec in CFG['metrics'].items():
    s=score_for_auc(c[spec['column']],spec['direction']); cc.append({'method':method,'n_pairs':len(c),'n_positive':int(c.y.sum()),'n_negative':int((1-c.y).sum()),'n_genera':c.genus_eval.nunique(),'roc_auc_micro':roc_auc_score(c.y,s),'average_precision':average_precision_score(c.y,s)})
pd.DataFrame(cc).to_csv(OUT/'common_complete_case_method_comparison.tsv',sep='\t',index=False)
# Mash sentinel sensitivity: exclude D==1 no-detectable-similarity outputs.
m=df[df.mash_distance.notna() & df.mash_distance.lt(1.0) & df.genus_eval.notna()].copy(); both=m.groupby('genus_eval').y.nunique(); m=m[m.genus_eval.isin(both[both==2].index)]
roc,ap,pr=roc_pr(m.y,m.mash_distance,'lower'); pred=threshold_pred(m.mash_distance,.05,'<='); mashrow={'analysis':'Mash D<1 detectable-only sensitivity','n_pairs':len(m),'n_positive':int(m.y.sum()),'n_negative':int((1-m.y).sum()),'n_genera':m.genus_eval.nunique(),'roc_auc_micro':roc,'average_precision':ap,'pr_auc_trapezoid':pr,**confusion_metrics(m.y,pred)}
pd.DataFrame([mashrow]).to_csv(OUT/'Mash_distance1_sensitivity_analysis.tsv',sep='\t',index=False)
print(pd.DataFrame(summary)[['method','n_pairs','n_positive','n_negative','n_genera','roc_auc_micro','average_precision','full_data_youden_cutoff_descriptive_only']].to_string(index=False))
print('\nCommon complete case\n',pd.DataFrame(cc).to_string(index=False))
print('\nMash sentinel sensitivity\n',pd.DataFrame([mashrow]).to_string(index=False))
