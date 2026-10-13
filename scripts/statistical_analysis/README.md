# Statistical analysis

Compact scripts for the manuscript species-discrimination analyses: continuous-score ROC/PR performance and fixed literature-supported species thresholds.

## Input

The large pairwise score table is not stored in GitHub. Both analysis scripts expect a local tab-separated file named:

`data/species_discrimination_scores.tsv`

Required header:

```text
method	reference	genus	class	score
```

Each row represents one method/pair/reference combination. `method` must be `16S rRNA`, `FastANI`, `Mash`, or `FastAAI`; `reference` is `NCBI` or `GTDB`; `class` is `SS` (same species) or `SGDS` (same genus/different species). Use percent units for 16S rRNA and FastANI, Mash distance on its native 0–1 scale, and raw FastAAI/Jaccard similarity. Rows without a usable numeric score should be omitted.

For each method/reference, only genera containing at least one SS and one SGDS comparison are retained, matching the manuscript analysis.

## Run

```bash
python scripts/statistical_analysis/roc_pr.py
python scripts/statistical_analysis/thresholds.py
python scripts/statistical_analysis/plot.py
```

`roc_pr.py` uses NCBI labels by default; use `--reference GTDB` for the GTDB sensitivity analysis.

## Outputs

Outputs are written to `outputs/statistical_analysis/` and include ROC/PR summaries, curve points, literature-threshold performance, and the combined ROC/PR figure.

Fixed criteria are 16S rRNA identity ≥98.65%, FastANI ≥95% ANI, and Mash distance ≤0.05. No universal literature-supported cutoff is imposed on raw FastAAI/Jaccard similarity. Mash is sign-reversed internally for ROC/PR calculations because lower distance indicates greater similarity.
