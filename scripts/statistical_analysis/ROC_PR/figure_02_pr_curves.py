from pathlib import Path
import pandas as pd, matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]; T=ROOT/'outputs/tables'; O=ROOT/'outputs/figures'; O.mkdir(parents=True,exist_ok=True)
s=pd.read_csv(T/'primary_roc_pr_threshold_summary_ncbi.tsv',sep='\t').set_index('method')
fig,ax=plt.subplots(figsize=(7,6))
for method in s.index:
 d=pd.read_csv(T/'curve_points'/f'{method}_pr.tsv',sep='\t'); ax.plot(d.recall,d.precision,label=f"{method} (AP={s.loc[method,'average_precision']:.3f})")
ax.set(xlabel='Recall',ylabel='Precision',title='Species discrimination precision–recall curves'); ax.legend(fontsize=8); fig.tight_layout()
for ext in ['png','pdf']: fig.savefig(O/f'Figure_02_PR_curves.{ext}',dpi=300)
