# Taxonomic Resolution of 16S rRNA, FastANI, Mash, and FastAAI across 30,495 Type-Strain Genomes

## Inputs

- `data/gca_accessions_30495.txt` — frozen GCA accession.version list
- `data/final_clean_taxonomy_30495_UPDATED.tsv` — frozen NCBI-derived taxonomy
- `data/ncbi_assembly_metadata_key_columns_30495.tsv` — assembly metadata
- `data/taxonomy_fill_source.tsv` — source of completed ranks
- `data/GTDB_R232_30495_taxonomy_status.tsv` — GTDB assignment status
- Genome FASTA assemblies and large pairwise score files are supplied separately.

## Analysis order

1. Retrieve NCBI TaxIDs and taxonomy: `scripts/taxonomy/fetch_ncbi_taxids.sh`.
2. Check GenID taxonomy and accession set: `scripts/taxonomy/audit_dictionary.py`.
3. Inspect missing upper-rank assignments: `scripts/taxonomy/taxonomy_fill_safe.py`.
4. Predict 16S genes: `scripts/16s_rRNA/run_rnammer.sh`.
5. Extract qualifying copies and calculate recovery: `scripts/16s_rRNA/process_rnammer.py`.
6. Align sequences and compute conservation: `scripts/16s_rRNA/conservation.py`.
7. Calculate 16S pairwise identities: `scripts/16s_rRNA/run_vsearch.sh`.
8. Generate FastANI, Mash and FastAAI scores: `docs/30K_FastANI_Mash_FastAAI_Commands.md`. For FastANI blocks and Mash matrix conversion, use `docs/score_generation.md`.
9. Assign scores to taxonomic groups: `scripts/taxonomy/classify_pairwise_ranks.py`.
10. Stage pairwise and GTDB inputs under `scripts/statistical_analysis/final_species_discrimination/inputs/`.
11. Run `prepare_inputs.py`, `primary_roc_pr.py`, `gtdb_sensitivity.py`, and robustness scripts.
12. Generate plots with scripts under `scripts/sankey_plots/` and `scripts/other_plots/`.

Commands and formats: [workflow](docs/workflow.md), [genome comparison](docs/30K_FastANI_Mash_FastAAI_Commands.md), [score-generation utilities](docs/score_generation.md).
