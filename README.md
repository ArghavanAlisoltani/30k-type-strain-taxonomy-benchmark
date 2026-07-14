## **Taxonomic Resolution of 16S rRNA, FastANI, Mash, and FastAAI across 30,495 Prokaryotic Type-Strain Genomes**

This study evaluates four complementary sequence-based approaches—16S rRNA identity,
FastANI, Mash, and FastAAI/Jaccard—across a frozen collection of prokaryotic
type-strain genomes. The repository provides the accession list, standardized taxonomy,
compact NCBI assembly metadata, provenance notes, and custom scripts.

## Dataset snapshot

The benchmark contains **30,495 unique GenBank assembly accessions**:

| Category | Count | Percentage |
|---|---:|---:|
| Bacteria | 29,383 | 96.4% |
| Archaea | 1,112 | 3.6% |
| Complete Genome | 5,970 | 19.6% |
| Chromosome | 480 | 1.6% |
| Scaffold | 9,676 | 31.7% |
| Contig | 14,369 | 47.1% |

Complete Genome and Chromosome assemblies together account for 6,450 genomes (21.2%).

## Data retrieval

Type-strain records were retrieved from the NCBI Datasets Genome web interface on
**4 May 2026** using the **type material only** filter:

- Bacteria: <https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2&type_material_only=true>
- Archaea: <https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2157&type_material_only=true>

GenBank assembly accessions beginning with `GCA_` were retained as primary identifiers.


## Main data files

- `data/final_clean_taxonomy_30495.csv`: standardized nine-rank taxonomy,
  one row per GCA accession.
- `data/ncbi_assembly_metadata_key_columns_30495.tsv`: compact 31-column subset
  of the original 150-column parsed NCBI metadata table.
- `data/gca_accessions_30495.txt`: one GCA accession per line.


## Citation



## Licensing


