# Taxonomic Resolution of 16S rRNA, FastANI, Mash, and FastAAI across 30,495 Type-Strain Genomes

## Inputs

- `data/gca_accessions_30495.txt` — frozen GCA accession.version list
- `data/final_clean_taxonomy_30495_UPDATED.tsv` — frozen NCBI-derived taxonomy
- `data/ncbi_assembly_metadata_key_columns_30495.tsv` — assembly metadata
- `data/taxonomy_fill_source.tsv` — source of completed ranks
- `data/GTDB_R232_30495_taxonomy_status.tsv` — GTDB assignment status
- Genome FASTA assemblies and large pairwise score files are supplied separately.

## Analysis order

1. Retrieve NCBI accession–TaxID mappings: `scripts/taxonomy/fetch_ncbi_taxids.sh`.
2. Inspect optional upper-rank fill proposals: `scripts/taxonomy/taxonomy_fill_safe.py` (dry-run default).
3. Predict bacterial and archaeal 16S rRNA genes: `scripts/16s_rRNA/run_rnammer.sh`.
4. Extract full-length copies, select one sequence per genome and summarize recovery: `scripts/16s_rRNA/process_rnammer.py`.
5. Align full-length 16S sequences with MAFFT and calculate site conservation: `scripts/16s_rRNA/conservation.py`.
6. Calculate 16S pairwise identity with VSEARCH: `scripts/16s_rRNA/run_vsearch.sh`.
7. Generate FastANI, Mash and FastAAI scores using [`docs/30K_FastANI_Mash_FastAAI_Commands.md`](docs/30K_FastANI_Mash_FastAAI_Commands.md).
8. Assign scores to taxonomic relationships: `scripts/taxonomy/classify_pairwise_ranks.py`.
9. Supply original pairwise and GTDB inputs under `scripts/statistical_analysis/final_species_discrimination/inputs/` (see `scripts/statistical_analysis/README.md`).
10. Run `prepare_inputs.py`, `primary_roc_pr.py`, `gtdb_sensitivity.py` and robustness scripts under `scripts/statistical_analysis/final_species_discrimination/scripts/`.
11. Run figure scripts under `scripts/sankey_plots/` and `scripts/other_plots/`.

Commands, file formats, and the two distinct VSEARCH modes are in [`docs/workflow.md`](docs/workflow.md). Code checks: `python -m unittest discover -s tests -v`.
