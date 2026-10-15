#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, sys, json, time
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'scripts'))
from common import pair_key_vec
RAW=ROOT/'inputs/raw_uploads'; POS=ROOT/'inputs/extracted/positive'; NEG=ROOT/'inputs/extracted/genus_negative'; DER=ROOT/'inputs/derived'; OUT=ROOT/'outputs/tables'; DER.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)

t0=time.time()
def stage(msg):
    print(f'[{time.time()-t0:8.1f}s] {msg}', flush=True)

stage('Loading frozen NCBI taxonomy')
tax=pd.read_csv(RAW/'final_clean_taxonomy_30495.csv',dtype=str,low_memory=False, usecols=['accession','phylum','genus','species'])
assert len(tax)==30495 and tax.accession.nunique()==30495, 'Frozen NCBI taxonomy must have 30,495 unique accessions'
tax=tax.set_index('accession',drop=False)
stage('Loading NCBI-GTDB R232 provenance master')
gtdb=pd.read_csv(RAW/'FINAL_30495_NCBI_GTDB_R232_provenance_master.tsv',sep='\t',dtype=str,low_memory=False, usecols=['GCA','GTDB_genus','GTDB_species','GTDB_status','genus_normalized_status','species_normalized_status','genus_exact_status','species_exact_status'])
assert len(gtdb)==30495 and gtdb.GCA.nunique()==30495, 'GTDB provenance master must have 30,495 unique GCA rows'
gtdb=gtdb.set_index('GCA',drop=False)
gtdb['_norm_consensus']=gtdb.genus_normalized_status.eq('AGREE') & gtdb.species_normalized_status.eq('AGREE')
gtdb['_strict_consensus']=gtdb.genus_exact_status.eq('AGREE') & gtdb.species_exact_status.eq('AGREE')

def add_ncbi_meta(df):
    # Only ranks used by this analysis are materialized; the full frozen taxonomy remains in inputs/raw_uploads.
    for side in [1,2]:
        gcol='genome'+str(side)
        for r in ['phylum','genus','species']:
            df[f'ncbi_{r}{side}']=df[gcol].map(tax[r])
    return df

def add_gtdb_meta(df):
    # Map only fields required for the fixed-pair GTDB sensitivity analysis.
    for side in [1,2]:
        gcol='genome'+str(side)
        df[f'GTDB_genus{side}']=df[gcol].map(gtdb['GTDB_genus'])
        df[f'GTDB_species{side}']=df[gcol].map(gtdb['GTDB_species'])
        df[f'name_consensus_normalized_genome{side}']=df[gcol].map(gtdb['_norm_consensus']).fillna(False).astype(bool)
        df[f'name_consensus_strict_genome{side}']=df[gcol].map(gtdb['_strict_consensus']).fillna(False).astype(bool)
    complete=df[['GTDB_genus1','GTDB_genus2','GTDB_species1','GTDB_species2']].notna().all(axis=1)
    sameg=df.GTDB_genus1.eq(df.GTDB_genus2); sames=df.GTDB_species1.eq(df.GTDB_species2)
    df['y_gtdb']=np.where(complete & sameg & sames,1,np.where(complete & sameg & ~sames,0,np.nan))
    df['gtdb_genus_eval']=np.where(pd.notna(df.y_gtdb),df.GTDB_genus1,np.nan)
    df['pair_relation_consensus']=pd.notna(df.y_gtdb) & df.y_ncbi.eq(df.y_gtdb)
    df['pair_name_consensus_normalized']=df.name_consensus_normalized_genome1 & df.name_consensus_normalized_genome2
    df['pair_name_consensus_strict']=df.name_consensus_strict_genome1 & df.name_consensus_strict_genome2
    return df

def canonicalize(df,a='genome1',b='genome2'):
    df=df.copy(); df['pair_id']=pair_key_vec(df[a],df[b]).values; return df

def assert_unique_subset(x,universe,name):
    if x.pair_id.duplicated().any(): raise AssertionError(f'{name}: duplicate canonical pairs in retained universe')
    extra=~x.pair_id.isin(universe)
    if extra.any(): raise AssertionError(f'{name}: {int(extra.sum())} retained rows outside intended pair universe')

stage('Preparing same-species pair universe')
# ---------- positive / same-species intended universe ----------
pos_mash=pd.read_csv(POS/'MASH.dat',sep='\t',header=None,names=['genome1','genome2','mash_distance'],dtype={'genome1':str,'genome2':str,'mash_distance':float})
pos_mash=canonicalize(pos_mash)
assert len(pos_mash)==16402 and pos_mash.pair_id.nunique()==16402
pos=pos_mash.copy(); pos=add_ncbi_meta(pos)
assert pos.ncbi_species1.eq(pos.ncbi_species2).all(), 'Positive Mash universe is not entirely same-species'
pos['y_ncbi']=1; pos['genus_eval']=pos.ncbi_genus1; pos['phylum_eval']=pos.ncbi_phylum1; pos['relationship_ncbi']='same_species'
pos_ids=set(pos.pair_id)

# FastANI positive
x=pd.read_csv(POS/'F-ANI.dat',sep='\t',header=None,names=['genome1','genome2','fastani'],dtype={'genome1':str,'genome2':str,'fastani':float}); x=canonicalize(x); assert_unique_subset(x,pos_ids,'positive FastANI'); pos=pos.merge(x[['pair_id','fastani']],on='pair_id',how='left')
# FastAAI positive
x=pd.read_csv(POS/'F-AAI.dat',sep='\t',dtype=str,keep_default_na=False).rename(columns={'#target':'genome1','query':'genome2','Jaccard':'fastaai_jaccard','SCPs':'fastaai_scps','ANI_Score':'fastaai_estimated_aai_raw'}); x=canonicalize(x); assert_unique_subset(x,pos_ids,'positive FastAAI'); x['fastaai_jaccard']=pd.to_numeric(x.fastaai_jaccard.replace({'N/A':np.nan,'NA':np.nan,'':np.nan}),errors='coerce'); x['fastaai_scps']=pd.to_numeric(x.fastaai_scps,errors='coerce'); pos=pos.merge(x[['pair_id','fastaai_jaccard','fastaai_scps','fastaai_estimated_aai_raw']],on='pair_id',how='left')
# 16S positive; quarantine rows outside intended pair universe
x=pd.read_csv(POS/'VSEARCH.dat',sep='\t',header=None,names=['genome1','genome2','identity_16s'],dtype={'genome1':str,'genome2':str,'identity_16s':float}); x=canonicalize(x)
extra=x[~x.pair_id.isin(pos_ids)].copy(); extra.to_csv(OUT/'QC_positive_16S_rows_outside_intended_universe.tsv',sep='\t',index=False)
x=x[x.pair_id.isin(pos_ids)].copy(); assert_unique_subset(x,pos_ids,'positive 16S'); pos=pos.merge(x[['pair_id','identity_16s']],on='pair_id',how='left')
pos=add_gtdb_meta(pos)
pos['mash_detectable_similarity']=pos.mash_distance.lt(1.0)

stage('Preparing same-genus/different-species pair universe')
# ---------- negative / same-genus-different-species intended universe ----------
neg_mash=pd.read_csv(NEG/'MASH.dat',sep='\t',header=None,names=['genome1','genome2','mash_distance'],dtype={'genome1':str,'genome2':str,'mash_distance':float}); neg_mash=canonicalize(neg_mash)
assert len(neg_mash)==1754127 and neg_mash.pair_id.nunique()==1754127
neg=neg_mash.copy(); neg=add_ncbi_meta(neg)
sgds=neg.ncbi_genus1.eq(neg.ncbi_genus2) & ~neg.ncbi_species1.eq(neg.ncbi_species2)
assert sgds.all(), 'Negative Mash universe is not entirely same-genus/different-species'
neg['y_ncbi']=0; neg['genus_eval']=neg.ncbi_genus1; neg['phylum_eval']=neg.ncbi_phylum1; neg['relationship_ncbi']='same_genus_diff_species'; neg_ids=set(neg.pair_id)
# FastANI negative with alignment fraction
x=pd.read_csv(NEG/'F-ANI.dat',sep='\t',header=None,names=['genome1','genome2','fastani','fastani_matched_fragments','fastani_total_fragments'],dtype={'genome1':str,'genome2':str,'fastani':float,'fastani_matched_fragments':float,'fastani_total_fragments':float}); x=canonicalize(x); assert_unique_subset(x,neg_ids,'negative FastANI'); x['fastani_alignment_fraction']=x.fastani_matched_fragments/x.fastani_total_fragments; neg=neg.merge(x[['pair_id','fastani','fastani_matched_fragments','fastani_total_fragments','fastani_alignment_fraction']],on='pair_id',how='left')
# FastAAI negative: no header, five columns
x=pd.read_csv(NEG/'F-AAI.dat',sep='\t',header=None,names=['genome1','genome2','fastaai_jaccard','fastaai_scps','fastaai_estimated_aai'],dtype=str,keep_default_na=False); x=canonicalize(x); assert_unique_subset(x,neg_ids,'negative FastAAI');
for c in ['fastaai_jaccard','fastaai_scps','fastaai_estimated_aai']: x[c]=pd.to_numeric(x[c].replace({'N/A':np.nan,'NA':np.nan,'':np.nan}),errors='coerce')
neg=neg.merge(x[['pair_id','fastaai_jaccard','fastaai_scps','fastaai_estimated_aai']],on='pair_id',how='left')
# 16S negative: original file contains a small noncanonical set; quarantine it and keep exact intended universe only
x=pd.read_csv(NEG/'16S.dat',sep='\t',header=None,names=['genome1','genome2','identity_16s'],dtype={'genome1':str,'genome2':str,'identity_16s':float}); x=canonicalize(x)
extra=x[~x.pair_id.isin(neg_ids)].copy(); extra=add_ncbi_meta(extra); extra.to_csv(OUT/'QC_negative_16S_rows_outside_intended_universe.tsv',sep='\t',index=False)
x=x[x.pair_id.isin(neg_ids)].copy(); assert_unique_subset(x,neg_ids,'negative 16S'); neg=neg.merge(x[['pair_id','identity_16s']],on='pair_id',how='left')
neg=add_gtdb_meta(neg); neg['mash_detectable_similarity']=neg.mash_distance.lt(1.0)

stage('Joining method scores and writing compact binary analysis table')
# ---------- save compact derived pair table ----------
# Raw uploads and extracted method files are retained separately, so the derived table keeps only
# fields required by the statistical pipeline rather than duplicating the full taxonomy master.
metriccols=['identity_16s','fastani','mash_distance','fastaai_jaccard','fastani_matched_fragments','fastani_total_fragments','fastani_alignment_fraction','fastaai_scps','fastaai_estimated_aai','fastaai_estimated_aai_raw','mash_detectable_similarity']
gtdbcols=['y_gtdb','gtdb_genus_eval','pair_relation_consensus','pair_name_consensus_normalized','pair_name_consensus_strict']
corecols=['pair_id','genome1','genome2','relationship_ncbi','y_ncbi','genus_eval','phylum_eval']+metriccols+gtdbcols
for c in corecols:
    if c not in pos: pos[c]=np.nan
    if c not in neg: neg[c]=np.nan
combined=pd.concat([pos[corecols],neg[corecols]],ignore_index=True)
combined.to_csv(DER/'binary_analysis_core.tsv',sep='\t',index=False)
stage('Derived binary_analysis_core.tsv written')

# ---------- audit summaries ----------
rows=[]
for label,z in [('same_species',pos),('same_genus_diff_species',neg)]:
    for method,col in [('16S','identity_16s'),('FastANI','fastani'),('Mash_numeric','mash_distance'),('Mash_detectable','mash_detectable_similarity'),('FastAAI_Jaccard','fastaai_jaccard')]:
        if method=='Mash_detectable': valid=z.mash_detectable_similarity.fillna(False).astype(bool)
        else: valid=z[col].notna()
        rows.append({'class':label,'method':method,'intended_pairs':len(z),'valid_or_detectable_pairs':int(valid.sum()),'unavailable_or_censored_pairs':int((~valid).sum()),'availability_pct':100*valid.mean()})
pd.DataFrame(rows).to_csv(OUT/'input_availability_by_class.tsv',sep='\t',index=False)

# GTDB master audit
pd.DataFrame([
    {'metric':'rows','value':len(gtdb)}, {'metric':'unique_GCA','value':gtdb.GCA.nunique()},
    {'metric':'GTDB_CLASSIFIED','value':int((gtdb.GTDB_status=='CLASSIFIED').sum())}, {'metric':'R232_QC_FAILED','value':int((gtdb.GTDB_status=='R232_QC_FAILED').sum())}, {'metric':'UNCLASSIFIED','value':int((gtdb.GTDB_status=='UNCLASSIFIED').sum())},
    {'metric':'genus_species_normalized_AGREE_genomes','value':int((gtdb.genus_normalized_status.eq('AGREE') & gtdb.species_normalized_status.eq('AGREE')).sum())},
    {'metric':'genus_species_strict_AGREE_genomes','value':int((gtdb.genus_exact_status.eq('AGREE') & gtdb.species_exact_status.eq('AGREE')).sum())}
]).to_csv(OUT/'GTDB_master_audit.tsv',sep='\t',index=False)
# pair label transition
eligible=combined[combined.y_gtdb.notna()].copy(); ct=pd.crosstab(eligible.y_ncbi,eligible.y_gtdb).rename_axis(index='NCBI_y',columns='GTDB_y').stack().reset_index(name='pair_count'); ct.to_csv(OUT/'NCBI_GTDB_pair_label_transition.tsv',sep='\t',index=False)
# FastAAI censoring audit
fa_pos=pos.fastaai_estimated_aai_raw.astype(str); censor=fa_pos.str.contains('>90',regex=False,na=False)
pd.DataFrame([{'class':'same_species','rows':len(pos),'continuous_estimated_AAI_available':int(pd.to_numeric(pos.fastaai_estimated_aai_raw,errors='coerce').notna().sum()),'censored_gt90_rows':int(censor.sum())}, {'class':'same_genus_diff_species','rows':len(neg),'continuous_estimated_AAI_available':int(neg.fastaai_estimated_aai.notna().sum()),'censored_gt90_rows':0}]).to_csv(OUT/'FastAAI_estimated_AAI_censoring_audit.tsv',sep='\t',index=False)
# overall input inventory
inventory=pd.DataFrame([
 {'file':'positive/MASH.dat','rows':len(pos_mash),'role':'same-species intended universe + Mash score'},
 {'file':'genus_negative/MASH.dat','rows':len(neg_mash),'role':'same-genus/different-species intended universe + Mash score'},
 {'file':'positive/F-ANI.dat','rows':pos.fastani.notna().sum(),'role':'same-species FastANI'},
 {'file':'genus_negative/F-ANI.dat','rows':neg.fastani.notna().sum(),'role':'genus-negative FastANI + fragment counts'},
 {'file':'positive/F-AAI.dat','rows':int(pos['fastaai_estimated_aai_raw'].notna().sum()),'role':'same-species FastAAI rows retained in intended universe; raw Jaccard; estimated AAI partly censored'},
 {'file':'genus_negative/F-AAI.dat','rows':int(neg['fastaai_scps'].notna().sum()),'role':'genus-negative FastAAI rows retained in intended universe; raw Jaccard + estimated AAI'},
 {'file':'positive/VSEARCH.dat','rows':int(pos['identity_16s'].notna().sum()+len(pd.read_csv(OUT/'QC_positive_16S_rows_outside_intended_universe.tsv',sep='\t'))),'role':'same-species 16S; noncanonical rows quarantined'},
 {'file':'genus_negative/16S.dat','rows':int(neg['identity_16s'].notna().sum()+len(extra)),'role':'genus-negative 16S; noncanonical rows quarantined'},
 {'file':'NCBI taxonomy','rows':len(tax),'role':'frozen reference taxonomy'},
 {'file':'NCBI-GTDB provenance master','rows':len(gtdb),'role':'GTDB R232 fixed-pair sensitivity mapping'}
]); inventory.to_csv(OUT/'input_file_inventory.tsv',sep='\t',index=False)
stage('Audit tables written')
print('Prepared:',len(pos),'positive and',len(neg),'negative intended pairs; combined',len(combined))
print(pd.DataFrame(rows).to_string(index=False))
