#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *

def auc_mask_presorted(score,y,mask):
    idx=np.argsort(score,kind='mergesort'); s=score[idx]; yy=y[idx]; mm=mask[idx]
    s=s[mm]; yy=yy[mm]
    if yy.size==0 or yy.min()==yy.max(): return np.nan
    starts=np.r_[0,np.flatnonzero(s[1:]!=s[:-1])+1]
    pos=np.add.reduceat((yy==1).astype(np.int64),starts); neg=np.add.reduceat((yy==0).astype(np.int64),starts)
    before=np.cumsum(neg)-neg; den=pos.sum()*neg.sum(); return float(np.sum(pos*before+0.5*pos*neg)/den)

p,n=load_method('FastAAI_Jaccard'); z=pd.concat([p,n],ignore_index=True).dropna(subset=['score','scps','genus_eval']).copy(); both=z.groupby('genus_eval',sort=False).y.nunique(); z=z[z.genus_eval.isin(both[both==2].index)].copy()
score=z.score.to_numpy(float); y=z.y.to_numpy(np.int8); scp=z.scps.to_numpy(float)
subsets=[('all_numeric_scp',np.ones(len(z),bool)),('SCP_ge_60',scp>=60),('SCP_ge_70',scp>=70),('SCP_ge_75',scp>=75),('SCP_ge_78',scp>=78),('SCP_lt_60',scp<60),('SCP_60_69',(scp>=60)&(scp<70)),('SCP_70_74',(scp>=70)&(scp<75)),('SCP_75_plus',scp>=75)]
rows=[]
for lab,sel in subsets:
 q=z.loc[sel]
 if q.y.nunique()<2: continue
 rows.append({'subset':lab,'n_pairs':len(q),'positive_pairs':int(q.y.sum()),'negative_pairs':int((q.y==0).sum()),'genera':q.genus_eval.nunique(),'median_SCPs':q.scps.median(),'roc_auc':auc_mask_presorted(score,y,sel)})
pd.DataFrame(rows).to_csv(OUT/'FastAAI_SCP_sensitivity.tsv',sep='\t',index=False)
bins=[0,59,69,74,77,79,91]; labels=['<=59','60-69','70-74','75-77','78-79','80-91']; z['scp_bin']=pd.cut(z.scps,bins=bins,labels=labels,include_lowest=True)
d=z.groupby(['y','scp_bin'],observed=True).agg(n=('score','size'),median_score=('score','median'),mean_score=('score','mean'),median_scps=('scps','median')).reset_index(); d.to_csv(OUT/'FastAAI_score_by_SCP_bin_and_class.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
