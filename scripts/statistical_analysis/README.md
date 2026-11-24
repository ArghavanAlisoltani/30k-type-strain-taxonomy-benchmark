# Species-discrimination analysis

Working directory: `scripts/statistical_analysis/final_species_discrimination/`.

## Required local inputs

- `inputs/raw_uploads/final_clean_taxonomy_30495.csv` — local CSV with `accession`, `phylum`, `genus`, `species` columns derived from frozen taxonomy
- `inputs/raw_uploads/FINAL_30495_NCBI_GTDB_R232_provenance_master.tsv` — local GTDB provenance master
- `inputs/extracted/positive/{MASH.dat,F-ANI.dat,F-AAI.dat,VSEARCH.dat}`
- `inputs/extracted/genus_negative/{MASH.dat,F-ANI.dat,F-AAI.dat,16S.dat}`

The original large numerical files and original full GTDB provenance master are not committed. The main repository taxonomy file is `data/final_clean_taxonomy_30495_UPDATED.tsv` (TSV), not the local CSV input expected by `prepare_inputs.py`.

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
