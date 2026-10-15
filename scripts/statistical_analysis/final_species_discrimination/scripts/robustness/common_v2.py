from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
EX = ROOT/'inputs'/'extracted'; RAW=ROOT/'inputs'/'raw_uploads'
OUT=ROOT/'outputs'/'robustness_tables'; FIG=ROOT/'outputs'/'robustness_figures'
OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
METHOD_INFO={'16S':('identity_16s',True),'FastANI':('fastani',True),'Mash':('mash_distance',False),'FastAAI_Jaccard':('fastaai_jaccard',True)}
_TAX=None; _MAPS=None

def taxonomy():
 global _TAX
 if _TAX is None: _TAX=pd.read_csv(RAW/'final_clean_taxonomy_30495.csv',dtype=str)
 return _TAX

def master(): return pd.read_csv(RAW/'FINAL_30495_NCBI_GTDB_R232_provenance_master.tsv',sep='\t',dtype=str)

def tax_maps():
 global _MAPS
 if _MAPS is None:
  t=taxonomy().set_index('accession'); _MAPS={c:t[c].to_dict() for c in ['phylum','order','family','genus','species']}
 return _MAPS

def canonicalize(df,a='g1',b='g2'):
 df=df.copy(); aa=df[a].astype(str); bb=df[b].astype(str); swap=aa.gt(bb)
 df['genome1']=aa.where(~swap,bb); df['genome2']=bb.where(~swap,aa); df['pair_id']=df.genome1+'||'+df.genome2
 return df

def add_tax(df,maps=None):
 if maps is None: maps=tax_maps()
 z=df.copy()
 for c in ['genus','species']:
  mp=maps[c]; z[c+'1']=z.genome1.map(mp); z[c+'2']=z.genome2.map(mp)
 z['genus_eval']=z.genus1
 # Negative/positive pairs are within genus/species, so genome1's higher taxonomy is the blocking label.
 z['family_eval']=z.genome1.map(maps['family']); z['order_eval']=z.genome1.map(maps['order']); z['phylum_eval']=z.genome1.map(maps['phylum'])
 return z

def load_method(method):
 maps=tax_maps(); valid=set(maps['species'])
 if method=='16S':
  p=pd.read_csv(EX/'positive'/'VSEARCH.dat',sep='\t',header=None,names=['g1','g2','score'],dtype={'g1':str,'g2':str})
  n=pd.read_csv(EX/'genus_negative'/'16S.dat',sep='\t',header=None,names=['g1','g2','score'],dtype={'g1':str,'g2':str})
 elif method=='FastANI':
  p=pd.read_csv(EX/'positive'/'F-ANI.dat',sep='\t',header=None,names=['g1','g2','score'],dtype={'g1':str,'g2':str})
  n=pd.read_csv(EX/'genus_negative'/'F-ANI.dat',sep='\t',header=None,names=['g1','g2','score','matched_fragments','total_fragments'],dtype={'g1':str,'g2':str})
  n['matched_fragments']=pd.to_numeric(n.matched_fragments,errors='coerce'); n['total_fragments']=pd.to_numeric(n.total_fragments,errors='coerce'); n['alignment_fraction']=n.matched_fragments/n.total_fragments
 elif method=='Mash':
  p=pd.read_csv(EX/'positive'/'MASH.dat',sep='\t',header=None,names=['g1','g2','score'],dtype={'g1':str,'g2':str})
  n=pd.read_csv(EX/'genus_negative'/'MASH.dat',sep='\t',header=None,names=['g1','g2','score'],dtype={'g1':str,'g2':str})
 elif method=='FastAAI_Jaccard':
  p=pd.read_csv(EX/'positive'/'F-AAI.dat',sep='\t',comment='#',header=None,names=['g1','g2','score','scps','estimated_aai'],dtype={'g1':str,'g2':str})
  n=pd.read_csv(EX/'genus_negative'/'F-AAI.dat',sep='\t',header=None,names=['g1','g2','score','scps','estimated_aai'],dtype={'g1':str,'g2':str})
  p['scps']=pd.to_numeric(p.scps,errors='coerce'); n['scps']=pd.to_numeric(n.scps,errors='coerce')
 else: raise ValueError(method)
 p=canonicalize(p); n=canonicalize(n); p['score']=pd.to_numeric(p.score,errors='coerce'); n['score']=pd.to_numeric(n.score,errors='coerce')
 p=p[p.genome1.isin(valid)&p.genome2.isin(valid)].copy(); n=n[n.genome1.isin(valid)&n.genome2.isin(valid)].copy()
 p=add_tax(p,maps); n=add_tax(n,maps)
 p=p[p.species1.eq(p.species2)].copy(); n=n[n.genus1.eq(n.genus2)&~n.species1.eq(n.species2)].copy(); p['y']=1; n['y']=0
 return p,n

def scored_primary(method):
 p,n=load_method(method); z=pd.concat([p,n],ignore_index=True).dropna(subset=['score','genus_eval']).copy(); cc=z.groupby('genus_eval',sort=False).y.nunique(); keep=set(cc[cc==2].index); return z[z.genus_eval.isin(keep)].copy()
def score_for_auc(z,method):
 s=z.score.to_numpy(float); return -s if method=='Mash' else s
def operating_rule(method,threshold,score): return score<=threshold if method=='Mash' else score>=threshold
