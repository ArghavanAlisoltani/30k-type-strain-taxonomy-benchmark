#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; D=ROOT/'outputs/tables/cv_by_method'; O=ROOT/'outputs/tables'
summary=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(D.glob('*_summary.tsv'))],ignore_index=True)
folds=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(D.glob('*_folds.tsv'))],ignore_index=True)
genus=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(D.glob('*_oof_genus_confusion.tsv'))],ignore_index=True)
summary.to_csv(O/'cv_summary_genus_blocked.tsv',sep='\t',index=False); folds.to_csv(O/'cv_fold_results_genus_blocked.tsv',sep='\t',index=False); genus.to_csv(O/'cv_oof_genus_confusion.tsv',sep='\t',index=False)
print(summary[['method','median_cv_empirical_cutoff','mean_heldout_roc_auc','sd_heldout_roc_auc','oof_sensitivity','oof_specificity']].to_string(index=False))
