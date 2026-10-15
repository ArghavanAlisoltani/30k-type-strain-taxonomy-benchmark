#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
M=master().set_index('GCA')['ncbi_assembly_level'].to_dict(); hi={'Complete Genome','Chromosome'}
def cat(a,b):
 x=M.get(a,'Unknown'); y=M.get(b,'Unknown')
 if x in hi and y in hi: return 'both_complete_or_chromosome'
 if x=='Contig' or y=='Contig': return 'any_contig'
 if x=='Scaffold' or y=='Scaffold': return 'any_scaffold_no_contig'
 return 'unknown_or_other'
def auc_selected(score,y,sel):
 s=score[sel]; yy=y[sel]; order=np.argsort(s,kind='mergesort');s=s[order];yy=yy[order]
 starts=np.r_[0,np.flatnonzero(s[1:]!=s[:-1])+1];pos=np.add.reduceat((yy==1).astype(np.int64),starts);neg=np.add.reduceat((yy==0).astype(np.int64),starts);before=np.cumsum(neg)-neg;return float(np.sum(pos*before+0.5*pos*neg)/(pos.sum()*neg.sum()))
rows=[]
for m in METHOD_INFO:
 z=scored_primary(m); z['assembly_category']=[cat(a,b) for a,b in zip(z.genome1,z.genome2)]; y=z.y.to_numpy(np.int8); score=score_for_auc(z,m)
 for c in ['both_complete_or_chromosome','any_scaffold_no_contig','any_contig','unknown_or_other']:
  sel=z.assembly_category.eq(c).to_numpy(); q=z.loc[sel]
  if len(q)==0 or q.y.nunique()<2: continue
  rows.append({'method':m,'assembly_category':c,'n_pairs':len(q),'positive_pairs':int(q.y.sum()),'negative_pairs':int((q.y==0).sum()),'genera':q.genus_eval.nunique(),'roc_auc':auc_selected(score,y,sel)})
pd.DataFrame(rows).to_csv(OUT/'assembly_level_roc_sensitivity.tsv',sep='\t',index=False); print(pd.DataFrame(rows).to_string(index=False))
