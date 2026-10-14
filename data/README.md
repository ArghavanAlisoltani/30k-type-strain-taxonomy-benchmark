# Data directory

Every accession-level file contains exactly **30,495 unique GCA accessions** or is derived from that fixed accession set.

- `final_clean_taxonomy_30495.csv`: standardized taxonomy.
- `taxonomy_fill_provenance.tsv`: rank-level record of repaired taxonomy fields; final values are synchronized to `final_clean_taxonomy_30495.csv`. `source` reports the provenance database recorded in the audit; `not_recorded` indicates no source was recorded for that repair in the audit workbook.
- `ncbi_assembly_metadata_key_columns_30495.tsv`: compact NCBI assembly metadata.
- `gca_accessions_30495.txt`: accession manifest.
- `dataset_summary.tsv`: domain and assembly-level summaries.

Large pairwise analysis tables are not stored in this repository because of their size. Required input formats are documented with the corresponding scripts under `scripts/`.
