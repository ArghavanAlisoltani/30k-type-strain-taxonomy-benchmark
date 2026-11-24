# Species-discrimination analysis

Working directory: `scripts/statistical_analysis/final_species_discrimination/`.

## Required local inputs

- `inputs/raw_uploads/final_clean_taxonomy_30495.csv` — local CSV with `accession`, `phylum`, `genus`, `species` columns derived from frozen taxonomy
- `inputs/raw_uploads/FINAL_30495_NCBI_GTDB_R232_provenance_master.tsv` — local GTDB provenance master
- `inputs/extracted/positive/{MASH.dat,F-ANI.dat,F-AAI.dat,VSEARCH.dat}`
- `inputs/extracted/genus_negative/{MASH.dat,F-ANI.dat,F-AAI.dat,16S.dat}`

Generate scores with [genome comparison commands](../../docs/30K_FastANI_Mash_FastAAI_Commands.md) or [score-generation utilities](../../docs/score_generation.md).

For `F-ANI.dat`, `prepare_inputs.py` reads `genome1`, `genome2`, `ANI` from `positive/` and `genome1`, `genome2`, `ANI`, `matched_fragments`, `total_query_fragments` from `genus_negative/`. Select the fixed SS and SGDS pair sets before staging these files.

Convert `data/final_clean_taxonomy_30495_UPDATED.tsv` to a CSV with `accession`, `phylum`, `genus`, `species` columns for `inputs/raw_uploads/final_clean_taxonomy_30495.csv`.

## Execution order

```bash
cd scripts/statistical_analysis/final_species_discrimination
python scripts/prepare_inputs.py
python scripts/primary_roc_pr.py
python scripts/gtdb_sensitivity.py
python scripts/genus_blocked_cv.py
python scripts/combine_genus_blocked_cv.py
```

Further scripts under `scripts/robustness/` use the same method-score inputs. Configuration: `config/analysis_config.json`. Core derived score table: `inputs/derived/binary_analysis_core.tsv`.
