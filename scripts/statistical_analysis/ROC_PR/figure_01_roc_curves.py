from pathlib import Path
import pandas as pd, matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]; T=ROOT/'outputs/tables'; O=ROOT/'outputs/figures'; O.mkdir(parents=True,exist_ok=True)
s=pd.read_csv(T/'primary_roc_pr_threshold_summary_ncbi.tsv',sep='\t').set_index('method')
fig,ax=plt.subplots(figsize=(7,6))
for method in s.index:
 d=pd.read_csv(T/'curve_points'/f'{method}_roc.tsv',sep='\t'); ax.plot(d.fpr,d.tpr,label=f"{method} (AUC={s.loc[method,'roc_auc_micro']:.4f})")
ax.plot([0,1],[0,1],linestyle='--',linewidth=1); ax.set(xlabel='False positive rate',ylabel='True positive rate',title='Species discrimination ROC: same species vs same genus/different species'); ax.legend(fontsize=8); fig.tight_layout()
for ext in ['png','pdf']: fig.savefig(O/f'Figure_01_ROC_curves.{ext}',dpi=300)
