# Statistical analysis

This directory contains the **exact analysis code used for the finalized species-discrimination results**. The canonical scripts are under `final_species_discrimination/`. No analysis data are stored in this repository.

The existing `roc_pr.py`, `thresholds.py`, and `plot.py` are compact convenience scripts created for easy application.

## Result-to-script map

| Result | Exact script(s) | Input type |
|---|---|---|
| Primary NCBI ROC/PR, literature-threshold performance, common complete-case analysis (Tables S4-S5) | `scripts/primary_roc_pr.py` | Derived pair-level table containing NCBI SS/SGDS labels and numeric 16S, FastANI, Mash, and FastAAI/Jaccard scores |
| GTDB R232 relabeling sensitivity (Table S4) | `scripts/gtdb_sensitivity.py` | Same derived pair-level table with GTDB genus/species labels and NCBI-GTDB concordance fields |
| Five-fold genus-blocked cross-validation (Table S7) | `scripts/genus_blocked_cv.py`, `scripts/combine_genus_blocked_cv.py` | Derived pair-level table; the combine script uses the per-method CV output tables |
| Species/species-pair collapsed analysis (Table S6A) | `scripts/robustness/species_pair_collapsed_roc.py` | Pairwise positive and congeneric-negative method outputs plus frozen NCBI taxonomy |
| Family- and order-blocked cross-validation (Table S7) | `scripts/robustness/family_order_blocked_cv.py`, `scripts/robustness/combine_family_order_blocked_cv.py` | Pairwise method outputs plus family/order taxonomic labels; the combine script uses per-method blocked-CV outputs |
| FastAAI single-copy-protein sensitivity (Table S9) | `scripts/robustness/fastaai_scp_sensitivity.py` | FastAAI positive and congeneric-negative pair outputs containing raw Jaccard score and recovered SCP count |
| Assembly-level sensitivity (Table S8) | `scripts/robustness/assembly_level_sensitivity.py` | Pairwise method outputs plus genome assembly-level metadata |
| Genus-cluster bootstrap AUC and paired complete-case AUC differences (Table S6B) | `scripts/robustness/genus_cluster_bootstrap_auc.py`, `scripts/robustness/paired_complete_case_bootstrap.py`, `scripts/robustness/combine_cluster_bootstrap.py` | Pairwise outputs for the four methods plus genus labels; the combine script uses the bootstrap output tables |

## Supporting files

- `config/analysis_config.json` stores the random seed, literature thresholds, score directions, and bootstrap/CV settings.
- `scripts/common.py` and `scripts/robustness/common_v2.py` provide shared analysis functions and input parsing.
- `scripts/prepare_inputs.py` builds the derived pair-level analysis table from the frozen taxonomy/provenance metadata and pairwise method outputs.

The large pairwise method outputs, taxonomy tables, and derived analysis tables are intentionally not stored in GitHub.

## Provenance

The Python files under `final_species_discrimination/` are verbatim copies from the canonical analysis bundle embedded in `MASTER_30K_ANALYSIS_AUDIT.ipynb`; only repository filenames/locations were changed for clarity.
