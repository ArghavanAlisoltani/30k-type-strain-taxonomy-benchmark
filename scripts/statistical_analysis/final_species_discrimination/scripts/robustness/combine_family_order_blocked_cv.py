#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common_v2 import *
fs=[]; ss=[]
for m in METHOD_INFO:
 fs.append(pd.read_csv(OUT/f'_blocked_cv_folds_{m}.tsv',sep='\t')); ss.append(pd.read_csv(OUT/f'_blocked_cv_summary_{m}.tsv',sep='\t'))
pd.concat(fs,ignore_index=True).to_csv(OUT/'phylogenetic_blocked_cv_folds_new_family_order.tsv',sep='\t',index=False)
pd.concat(ss,ignore_index=True).to_csv(OUT/'phylogenetic_blocked_cv_summary_new_family_order.tsv',sep='\t',index=False)
print(pd.concat(ss,ignore_index=True).to_string(index=False))
