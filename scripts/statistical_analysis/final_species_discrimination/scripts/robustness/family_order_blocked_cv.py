#!/usr/bin/env python3
import argparse,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
from sklearn.metrics import roc_auc_score, confusion_matrix
SEED=20260915

def assign_groups_balanced(df,gcol,k=5):
 gs=df.groupby(gcol,sort=False).y.agg(['size','sum']).reset_index().rename(columns={'sum':'pos'})
 rng=np.random.default_rng(SEED); gs['jitter']=rng.random(len(gs))
 posg=gs[gs.pos>0].sort_values(['pos','size','jitter'],ascending=[False,False,True]).reset_index(drop=True)
 negg=gs[gs.pos==0].sort_values(['size','jitter'],ascending=[False,True]).reset_index(drop=True)
 fp=np.zeros(k,float); fn=np.zeros(k,float); mapping={}
 for _,r in posg.iterrows():
  f=int(np.lexsort((np.arange(k),fn,fp))[0]); mapping[r[gcol]]=f; fp[f]+=r.pos; fn[f]+=r['size']
 for _,r in negg.iterrows():
  f=int(np.argmin(fn)); mapping[r[gcol]]=f; fn[f]+=r['size']
 return df[gcol].map(mapping).to_numpy(np.int8)

def exact_training_youden(decision,y,trainmask):
 order=np.argsort(-decision,kind='mergesort'); s=decision[order]; yy=y[order]; tm=trainmask[order]
 pos=(tm&(yy==1)).astype(np.int64); neg=(tm&(yy==0)).astype(np.int64); P=pos.sum(); N=neg.sum(); tp=np.cumsum(pos); fp=np.cumsum(neg)
 ends=np.r_[np.flatnonzero(s[1:]!=s[:-1]),len(s)-1]; J=tp[ends]/P-fp[ends]/N; return float(s[ends[int(np.argmax(J))]])

ap=argparse.ArgumentParser(); ap.add_argument('method',choices=list(METHOD_INFO)); args=ap.parse_args(); m=args.method
z=scored_primary(m); summ=[]; folds=[]
for level,gcol in [('family','family_eval'),('order','order_eval')]:
 q=z.dropna(subset=[gcol]).reset_index(drop=True).copy(); decision=(-q.score.to_numpy(float) if m=='Mash' else q.score.to_numpy(float)); y=q.y.to_numpy(np.int8); foldid=assign_groups_balanced(q,gcol,5); oof=np.zeros(len(q),bool); local=[]
 for f in range(5):
  test=foldid==f; train=~test; td=exact_training_youden(decision,y,train); native=-td if m=='Mash' else td; pred=decision[test]>=td; oof[test]=pred; yt=y[test]
  if np.unique(yt).size<2: raise RuntimeError(f'{m} {level} fold {f+1} lacks a class')
  tn,fp,fn,tp=confusion_matrix(yt,pred,labels=[0,1]).ravel(); auc=roc_auc_score(yt,decision[test]); rec={'method':m,'block_level':level,'fold':f+1,'threshold_from_training':native,'n_test':int(test.sum()),'n_test_groups':q.loc[test,gcol].nunique(),'roc_auc':auc,'TP':tp,'FN':fn,'TN':tn,'FP':fp,'sensitivity':tp/(tp+fn),'specificity':tn/(tn+fp)}; folds.append(rec); local.append(rec)
 ff=pd.DataFrame(local); tn,fp,fn,tp=confusion_matrix(y,oof,labels=[0,1]).ravel(); summ.append({'method':m,'block_level':level,'folds':5,'n_pairs':len(q),'n_groups':q[gcol].nunique(),'fold_assignment':'whole-group deterministic greedy balance of positive and total pair counts','mean_heldout_roc_auc':ff.roc_auc.mean(),'sd_heldout_roc_auc':ff.roc_auc.std(ddof=1),'min_heldout_roc_auc':ff.roc_auc.min(),'max_heldout_roc_auc':ff.roc_auc.max(),'median_training_youden_cutoff':ff.threshold_from_training.median(),'oof_sensitivity':tp/(tp+fn),'oof_specificity':tn/(tn+fp),'oof_balanced_accuracy':0.5*(tp/(tp+fn)+tn/(tn+fp))})
pd.DataFrame(folds).to_csv(OUT/f'_blocked_cv_folds_{m}.tsv',sep='\t',index=False); pd.DataFrame(summ).to_csv(OUT/f'_blocked_cv_summary_{m}.tsv',sep='\t',index=False); print(pd.DataFrame(summ).to_string(index=False))
