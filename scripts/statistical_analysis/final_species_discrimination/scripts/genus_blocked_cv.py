#!/usr/bin/env python3
from pathlib import Path
import argparse, pandas as pd, numpy as np, sys
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'scripts'))
from common import load_config, score_for_auc, youden_threshold, threshold_pred, confusion_metrics, genus_confusion_table, bootstrap_confusion_by_group
CFG=load_config(ROOT); parser=argparse.ArgumentParser(); parser.add_argument('--method',required=True); args=parser.parse_args(); method=args.method
if method not in CFG['metrics']: raise SystemExit(f'Unknown method {method}')
spec=CFG['metrics'][method]; col=spec['column']; IN=ROOT/'inputs/derived/binary_analysis_core.tsv'; OUT=ROOT/'outputs/tables/cv_by_method'; OUT.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(IN,sep='\t',usecols=['pair_id','y_ncbi','genus_eval',col],low_memory=False).rename(columns={'y_ncbi':'y'})
d=d[d[col].notna() & d.genus_eval.notna()].copy(); both=d.groupby('genus_eval').y.nunique(); d=d[d.genus_eval.isin(both[both==2].index)].reset_index(drop=True)
cv=StratifiedGroupKFold(n_splits=CFG['cv_folds'],shuffle=True,random_state=CFG['random_seed']); foldrows=[]; genusrows=[]
for fold,(tr,te) in enumerate(cv.split(d,d.y,groups=d.genus_eval),1):
    train=d.iloc[tr]; test=d.iloc[te]; thr,j=youden_threshold(train.y,train[col],spec['direction']); rule='>=' if spec['direction']=='higher' else '<='
    pred=threshold_pred(test[col],thr,rule); met=confusion_metrics(test.y,pred); aucv=float(roc_auc_score(test.y,score_for_auc(test[col],spec['direction']))); ap=float(average_precision_score(test.y,score_for_auc(test[col],spec['direction'])))
    foldrows.append({'method':method,'fold':fold,'train_genera':train.genus_eval.nunique(),'test_genera':test.genus_eval.nunique(),'train_pairs':len(train),'test_pairs':len(test),'train_positive':int(train.y.sum()),'test_positive':int(test.y.sum()),'empirical_youden_threshold':thr,'train_youden_J':j,'heldout_roc_auc':aucv,'heldout_average_precision':ap,**met})
    tmp=test[['genus_eval','y']].copy(); tmp['pred']=pred; gt=genus_confusion_table(tmp,'pred','genus_eval','y'); gt['method']=method; gt['fold']=fold; gt['threshold']=thr; genusrows.append(gt)
fr=pd.DataFrame(foldrows); gg=pd.concat(genusrows,ignore_index=True); totals=gg[['tp','fn','tn','fp']].sum(); tp,fn,tn,fp=[int(totals[k]) for k in ['tp','fn','tn','fp']]
# derive pooled metrics directly from totals
sens=tp/(tp+fn); specv=tn/(tn+fp); prec=tp/(tp+fp); npv=tn/(tn+fn); bal=(sens+specv)/2; f1=2*prec*sens/(prec+sens)
ci=bootstrap_confusion_by_group(gg.groupby('genus',as_index=False)[['tp','fn','tn','fp']].sum(),B=CFG['binary_cluster_bootstrap_replicates'],seed=CFG['random_seed'])
sumrow={'method':method,'n_folds':len(fr),'n_oof_pairs':int(fr.test_pairs.sum()),'n_unique_genera':int(gg.genus.nunique()),'median_cv_empirical_cutoff':float(fr.empirical_youden_threshold.median()),'cutoff_min':float(fr.empirical_youden_threshold.min()),'cutoff_max':float(fr.empirical_youden_threshold.max()),'mean_heldout_roc_auc':float(fr.heldout_roc_auc.mean()),'sd_heldout_roc_auc':float(fr.heldout_roc_auc.std(ddof=1)),'min_heldout_roc_auc':float(fr.heldout_roc_auc.min()),'max_heldout_roc_auc':float(fr.heldout_roc_auc.max()),'mean_heldout_average_precision':float(fr.heldout_average_precision.mean()),'oof_tp':tp,'oof_fn':fn,'oof_tn':tn,'oof_fp':fp,'oof_sensitivity':sens,'oof_specificity':specv,'oof_precision':prec,'oof_npv':npv,'oof_balanced_accuracy':bal,'oof_f1':f1,**{f'oof_{k}':v for k,v in ci.items()}}
fr.to_csv(OUT/f'{method}_folds.tsv',sep='\t',index=False); gg.to_csv(OUT/f'{method}_oof_genus_confusion.tsv',sep='\t',index=False); pd.DataFrame([sumrow]).to_csv(OUT/f'{method}_summary.tsv',sep='\t',index=False)
print(pd.DataFrame([sumrow]).to_string(index=False))
