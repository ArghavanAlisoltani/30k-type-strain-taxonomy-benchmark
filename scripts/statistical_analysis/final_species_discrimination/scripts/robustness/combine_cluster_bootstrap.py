#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
ss=[];rr=[]
for m in METHOD_INFO:
 ss.append(pd.read_csv(OUT/f'_cluster_bootstrap_summary_{m}.tsv',sep='\t'));d=pd.read_csv(OUT/f'_cluster_bootstrap_{m}.tsv',sep='\t');d.insert(0,'method',m);rr.append(d)
ss.append(pd.read_csv(OUT/'_paired_complete_case_bootstrap_summary.tsv',sep='\t'))
pd.concat(ss,ignore_index=True).to_csv(OUT/'cluster_bootstrap_auc_summary.tsv',sep='\t',index=False);pd.concat(rr,ignore_index=True).to_csv(OUT/'cluster_bootstrap_method_auc_replicates.tsv.gz',sep='\t',index=False,compression='gzip')
print(pd.concat(ss,ignore_index=True).to_string(index=False))
