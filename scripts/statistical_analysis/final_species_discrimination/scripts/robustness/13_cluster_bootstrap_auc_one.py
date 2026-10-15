#!/usr/bin/env python3
import argparse,os,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
from sklearn.metrics import roc_auc_score
N_BOOT=int(os.environ.get('N_BOOT','200')); SEED=int(os.environ.get('BOOT_SEED','20260915'))

def prepare(y,s,group):
 groups=pd.Index(sorted(pd.unique(group))); gi=pd.Categorical(group,categories=groups).codes
 order=np.argsort(s,kind='mergesort'); ss=s[order]; yy=y[order]; gg=gi[order]
 starts=np.r_[0,np.flatnonzero(ss[1:]!=ss[:-1])+1]
 return groups,yy,gg,starts

def weighted_auc(yy,gg,starts,mult):
 w=mult[gg].astype(np.float64); posw=w*(yy==1); negw=w*(yy==0)
 pg=np.add.reduceat(posw,starts); ng=np.add.reduceat(negw,starts); before=np.cumsum(ng)-ng; den=pg.sum()*ng.sum()
 return float(np.sum(pg*before+0.5*pg*ng)/den)

ap=argparse.ArgumentParser(); ap.add_argument('method',choices=list(METHOD_INFO)); args=ap.parse_args();m=args.method
z=scored_primary(m).dropna(subset=['score','genus_eval']).copy(); y=z.y.to_numpy(np.int8); s=score_for_auc(z,m); groups,yy,gg,starts=prepare(y,s,z.genus_eval.astype(str).to_numpy())
rng=np.random.default_rng(SEED+list(METHOD_INFO).index(m)); B=[]
for b in range(N_BOOT):
 mult=rng.multinomial(len(groups),np.repeat(1/len(groups),len(groups))); B.append(weighted_auc(yy,gg,starts,mult))
point=roc_auc_score(y,s); lo,hi=np.percentile(B,[2.5,97.5])
pd.DataFrame({'replicate':np.arange(1,N_BOOT+1),'roc_auc':B}).to_csv(OUT/f'_cluster_bootstrap_{m}.tsv',sep='\t',index=False)
pd.DataFrame([{'analysis':'genus_cluster_bootstrap_micro_auc','method_or_contrast':m,'point_estimate':point,'ci_low':lo,'ci_high':hi,'bootstrap_replicates':N_BOOT,'n_pairs':len(z),'n_genera':len(groups),'note':'Whole genera sampled with replacement; all pairs in a sampled genus receive the genus multiplicity as weight.'}]).to_csv(OUT/f'_cluster_bootstrap_summary_{m}.tsv',sep='\t',index=False)
print(m,point,lo,hi)
