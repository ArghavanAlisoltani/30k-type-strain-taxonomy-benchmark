# Statistical analysis

Compact scripts for the manuscript species-discrimination analyses: continuous-score ROC/PR performance and fixed literature-supported species thresholds.

## Input

Both analysis scripts use one tab-separated long-format table:

`data/species_discrimination_scores.tsv`

Required columns:

| column | meaning |
|---|---|
| `method` | `16S rRNA`, `FastANI`, `Mash`, or `FastAAI` |
| `reference` | taxonomy used for the labels, e.g. `NCBI` or `GTDB` |
| `genus` | genus under that reference taxonomy |
| `class` | `SS` (same species) or `SGDS` (same genus/different species) |
| `score` | numeric method score |

Use percent units for 16S rRNA identity and FastANI (e.g. `98.65`, `95.0`), Mash distance on its native 0–1 scale, and raw FastAAI/Jaccard similarity. Rows without a usable numeric score should be omitted. For each method/reference, the scripts retain only genera containing at least one SS and one SGDS comparison, matching the manuscript discrimination analysis.

## Run

From the repository root:

```bash
python scripts/statistical_analysis/roc_pr.py
python scripts/statistical_analysis/thresholds.py
python scripts/statistical_analysis/plot.py
```

`roc_pr.py` uses NCBI labels by default; another label set can be selected with `--reference GTDB`.

## Outputs

Files are written to `outputs/statistical_analysis/`:

- `roc_pr_summary.tsv` — eligible pair counts, ROC AUC, and average precision.
- `curve_points/*_roc.tsv` and `*_pr.tsv` — points used to draw the curves.
- `literature_thresholds.tsv` — TP/FN/TN/FP and sensitivity, specificity, precision, F1, and balanced accuracy.
- `roc_pr_curves.pdf` and `.png` — combined ROC and PR figure.

The fixed criteria are 16S rRNA identity ≥98.65%, FastANI ≥95% ANI, and Mash distance ≤0.05. No universal literature-supported cutoff is imposed on raw FastAAI/Jaccard similarity. For ROC/PR calculations, Mash is sign-reversed internally because lower distance indicates greater similarity; the other methods use their native score direction.
