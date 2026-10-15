#!/usr/bin/env python3
import os,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
from sklearn.metrics import roc_auc_score
N_BOOT=int(os.environ.get('N_BOOT','200')); SEED=int(os.environ.get('BOOT_SEED','20260915'))

def prep(y,s,group):
 groups=pd.Index(sorted(pd.unique(group)));gi=pd.Categorical(group,categories=groups).codes;order=np.argsort(s,kind='mergesort');ss=s[order];yy=y[order];gg=gi[order];starts=np.r_[0,np.flatnonzero(ss[1:]!=ss[:-1])+1];return groups,yy,gg,starts

def wauc(yy,gg,starts,mult):
 w=mult[gg].astype(float);pg=np.add.reduceat(w*(yy==1),starts);ng=np.add.reduceat(w*(yy==0),starts);before=np.cumsum(ng)-ng;return float(np.sum(pg*before+0.5*pg*ng)/(pg.sum()*ng.sum()))
# Build common pair set
frames=[]
for m in METHOD_INFO:
 p,n=load_method(m);q=pd.concat([p,n],ignore_index=True)[['pair_id','genus_eval','y','score']].dropna(subset=['score']).rename(columns={'score':m});frames.append(q)
base=frames[0]
for q in frames[1:]: base=base.merge(q[['pair_id',q.columns[-1]]],on='pair_id',how='inner')
both=base.groupby('genus_eval',sort=False).y.nunique();base=base[base.genus_eval.isin(both[both==2].index)].copy(); y=base.y.to_numpy(np.int8); group=base.genus_eval.astype(str).to_numpy()
prepared={}; point={}
for m in METHOD_INFO:
 s=-base[m].to_numpy(float) if m=='Mash' else base[m].to_numpy(float);prepared[m]=prep(y,s,group);point[m]=roc_auc_score(y,s)
groups=prepared['16S'][0]
# all prepared share same groups by same complete-case rows
rng=np.random.default_rng(SEED+999); contrasts=[('FastANI','FastAAI_Jaccard'),('Mash','FastAAI_Jaccard'),('FastANI','Mash'),('FastANI','16S'),('Mash','16S')]; vals={c:[] for c in contrasts}
for b in range(N_BOOT):
 mult=rng.multinomial(len(groups),np.repeat(1/len(groups),len(groups))); auc={}
 for m in METHOD_INFO:
  _,yy,gg,starts=prepared[m];auc[m]=wauc(yy,gg,starts,mult)
 for c in contrasts: vals[c].append(auc[c[0]]-auc[c[1]])
rows=[]; reps=[]
for a,c in contrasts:
 lo,hi=np.percentile(vals[(a,c)],[2.5,97.5]);rows.append({'analysis':'paired_complete_case_genus_cluster_bootstrap_auc_difference','method_or_contrast':f'{a}-{c}','point_estimate':point[a]-point[c],'ci_low':lo,'ci_high':hi,'bootstrap_replicates':N_BOOT,'n_pairs':len(base),'n_genera':len(groups),'note':'Identical complete-case pair set for all four methods; whole genera resampled.'});reps += [{'contrast':f'{a}-{c}','replicate':i+1,'auc_difference':v} for i,v in enumerate(vals[(a,c)])]
pd.DataFrame(rows).to_csv(OUT/'_paired_complete_case_bootstrap_summary.tsv',sep='\t',index=False);pd.DataFrame(reps).to_csv(OUT/'paired_complete_case_bootstrap_replicates.tsv.gz',sep='\t',index=False,compression='gzip');print(pd.DataFrame(rows).to_string(index=False))
