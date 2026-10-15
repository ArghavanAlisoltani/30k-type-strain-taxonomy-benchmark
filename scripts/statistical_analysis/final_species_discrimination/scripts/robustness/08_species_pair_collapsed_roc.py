#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
from sklearn.metrics import roc_auc_score, average_precision_score

rows=[]; details=[]
for m in METHOD_INFO:
    p,n=load_method(m)
    p=p.dropna(subset=['score','species1','genus1']).copy(); n=n.dropna(subset=['score','species1','species2','genus1']).copy()
    # one value per species for within-species class
    pc=(p.groupby(['genus1','species1'],as_index=False).score.median()
          .rename(columns={'genus1':'genus','species1':'unit1'}))
    pc['unit2']=pc.unit1; pc['y']=1
    # one value per unordered interspecies pair within genus
    spairs=np.sort(n[['species1','species2']].astype(str).to_numpy(),axis=1)
    n=n.copy(); n['unit1']=spairs[:,0]; n['unit2']=spairs[:,1]
    nc=(n.groupby(['genus1','unit1','unit2'],as_index=False).score.median()
          .rename(columns={'genus1':'genus'})); nc['y']=0
    c=pd.concat([pc[['genus','unit1','unit2','score','y']],nc[['genus','unit1','unit2','score','y']]],ignore_index=True)
    both=c.groupby('genus').y.nunique(); c=c[c.genus.isin(both[both==2].index)].copy()
    auc=roc_auc_score(c.y, -c.score if m=='Mash' else c.score)
    ap=average_precision_score(c.y, -c.score if m=='Mash' else c.score)
    rows.append({'method':m,'collapsed_units':len(c),'positive_species_units':int(c.y.sum()),'negative_species_pair_units':int((c.y==0).sum()),'eligible_genera':c.genus.nunique(),'roc_auc':auc,'average_precision':ap,'aggregation':'median score per species (positive) and per unordered species pair within genus (negative)'})
    c.insert(0,'method',m); details.append(c)
pd.DataFrame(rows).to_csv(OUT/'species_pair_collapsed_roc_summary.tsv',sep='\t',index=False)
pd.concat(details,ignore_index=True).to_csv(OUT/'species_pair_collapsed_units.tsv.gz',sep='\t',index=False,compression='gzip')
print(pd.DataFrame(rows).to_string(index=False))
