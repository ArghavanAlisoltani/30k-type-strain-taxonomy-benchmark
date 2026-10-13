# GTDB R232 analysis

This folder documents the final GTDB R232 workflow used in the revised manuscript.

## Workflow

1. The 30,495 manuscript GCA accession.version identifiers were matched exactly against the GTDB R232 bacterial and archaeal metadata using `ncbi_genbank_assembly_accession`.
2. **29,250** genomes mapped directly to GTDB R232 metadata; **1,245** did not.
3. Of the 1,245 unmatched genomes, **90** were listed as GTDB R232 QC failures and were retained as a separate provenance category.
4. The remaining **1,155** manuscript assemblies were submitted to GTDB-Tk v2.7.2 with the R232 reference package. **1,153** produced usable GTDB assignments; **2** could not be classified because of insufficient marker-sequence alignment.
5. The final rank-specific NCBI–GTDB comparison used only genomes with non-empty assignments in both taxonomies. At the domain level this yielded **30,403 comparable genomes**; the 92 non-comparable genomes were the 90 GTDB QC failures plus the 2 insufficient-MSA genomes.
6. Exact taxon-name agreement was the primary comparison. Removal of GTDB alphabetic placeholder suffixes (for example `_A` or `_B`) was evaluated only as a secondary sensitivity analysis.

## Required inputs

Repository inputs:

- `data/gca_accessions_30495.txt` — one exact GCA accession.version per line.
- `data/final_clean_taxonomy_30495.csv` — frozen NCBI-derived reference taxonomy.

External GTDB R232 inputs:

- `bac120_metadata_r232.tsv.gz`
- `ar53_metadata_r232.tsv.gz`
- `qc_failed_r232.tsv`
- GTDB-Tk v2.7.2 R232 reference package for assemblies absent from the metadata mapping.

Large GTDB intermediate/output tables are not stored in GitHub. A local merged table sufficient for the concordance analysis should contain at least:

```text
GCA	NCBI_domain	NCBI_phylum	NCBI_class	NCBI_order	NCBI_family	NCBI_genus	NCBI_species	GTDB_domain	GTDB_phylum	GTDB_class	GTDB_order	GTDB_family	GTDB_genus	GTDB_species	GTDB_source
```

`GTDB_source` should distinguish direct R232 metadata mapping, GTDB-Tk classification, R232 QC failure, and insufficient-MSA cases.

Standard GTDB-Tk commands/options are described in the manuscript Methods; this repository focuses on the manuscript-specific mapping, relabeling, and comparison logic.
