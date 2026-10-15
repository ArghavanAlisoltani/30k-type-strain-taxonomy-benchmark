# Statistical analysis

This directory contains the **exact analysis code used for the finalized species-discrimination results**. The canonical scripts are under `final_species_discrimination/`, which mirrors the original analysis-bundle layout so the script contents can remain unchanged. No analysis data are stored in this repository.

The existing `roc_pr.py`, `thresholds.py`, and `plot.py` are compact convenience scripts created later; they are not the canonical scripts used to generate the final reported values.

## Result-to-script map

| Result | Exact script(s) | Input type |
|---|---|---|
| Primary NCBI ROC/PR, literature-threshold performance, common complete-case analysis (Tables S4-S5) | `scripts/01_primary_roc_pr.py` | Derived pair-level table containing NCBI SS/SGDS labels and numeric 16S, FastANI, Mash, and FastAAI/Jaccard scores |
| GTDB R232 relabeling sensitivity (Table S4) | `scripts/04_gtdb_sensitivity.py` | Same derived pair-level table with GTDB genus/species labels and NCBI-GTDB concordance fields |
| Five-fold genus-blocked cross-validation (Table S7) | `scripts/02_genus_blocked_cv_one.py`, `scripts/02b_combine_cv_results.py` | Derived pair-level table; the combine script uses the per-method CV output tables |
| Species/species-pair collapsed analysis (Table S6A) | `scripts/robustness/08_species_pair_collapsed_roc.py` | Pairwise positive and congeneric-negative method outputs plus frozen NCBI taxonomy |
| Family- and order-blocked analyses (Table S7) | `scripts/robustness/09_phylogenetic_blocked_cv_one.py`, `scripts/robustness/09b_combine_blocked_cv.py` | Pairwise method outputs plus taxonomic lineage labels; the combine script uses per-method blocked-CV outputs |
| FastAAI single-copy-protein sensitivity (Table S9) | `scripts/robustness/11_fastaai_scp_sensitivity.py` | FastAAI positive and congeneric-negative pair outputs containing raw Jaccard score and recovered SCP count |
| Assembly-level sensitivity (Table S8) | `scripts/robustness/12_assembly_level_sensitivity.py` | Pairwise method outputs plus genome assembly-level metadata |
| Genus-cluster bootstrap AUC and paired complete-case AUC differences (Table S6B) | `scripts/robustness/13_cluster_bootstrap_auc_one.py`, `scripts/robustness/13b_paired_complete_case_bootstrap.py`, `scripts/robustness/13c_combine_cluster_bootstrap.py` | Pairwise outputs for the four methods plus genus labels; the combine script uses the bootstrap output tables |

## Supporting files

- `config/analysis_config.json` stores the random seed, literature thresholds, score directions, and bootstrap/CV settings.
- `scripts/common.py` and `scripts/robustness/common_v2.py` provide shared analysis functions and input parsing.
- `scripts/00_prepare_inputs.py` builds the derived pair-level analysis table from the frozen taxonomy/provenance metadata and pairwise method outputs.

The large pairwise method outputs, taxonomy tables, and derived analysis tables are intentionally not stored in GitHub.

## Provenance

The Python files under `final_species_discrimination/` are verbatim copies from the canonical analysis bundle embedded in `MASTER_30K_ANALYSIS_AUDIT.ipynb`; only their repository location is different. Numeric prefixes are retained to preserve direct correspondence with the original analysis workflow.
