# Workflow commands and input/output formats

All paths below are examples. Preserve the original 4 May 2026 accession manifest and frozen taxonomic labels. Running NCBI Datasets again retrieves current records, not necessarily the historical classifications.

## 1. NCBI TaxIDs and taxonomy

```bash
bash scripts/taxonomy/fetch_ncbi_taxids.sh /absolute/path/to/data/gca_accessions_30495.txt /absolute/path/to/new_taxonomy_output
```

Outputs: `assembly_summary.jsonl`, `assembly_taxids.tsv`, `taxids.txt`, `taxonomy.jsonl`. The final standardized taxonomy in `data/final_clean_taxonomy_30495_UPDATED.tsv` remains the benchmark reference.

`taxonomy_fill_safe.py` proposes missing kingdom–family ranks in `GenID={...}` files. It does not update genus/species or perform external LPSN/GTDB curation.

```bash
python scripts/taxonomy/taxonomy_fill_safe.py --input /absolute/path/to/original_GenID.txt
python scripts/taxonomy/taxonomy_fill_safe.py --input /absolute/path/to/original_GenID.txt --apply --output /absolute/path/to/new_GenID.txt --audit /absolute/path/to/proposed_fill.tsv
```

`classify_pairwise_ranks.py` accepts a taxonomy TSV with `accession` and `domain,kingdom,phylum,class,order,family,genus,species` columns, plus a score TSV with `genome1,genome2,score`. Rows without sufficient lower-rank information remain unresolved.

```bash
python scripts/taxonomy/classify_pairwise_ranks.py --taxonomy /absolute/path/to/frozen_taxonomy.tsv --scores /absolute/path/to/pairs.tsv --output /absolute/path/to/new_labeled_pairs.tsv
```

## 2. RNAmmer v1.2 (Bacteria and Archaea)

```bash
bash scripts/16s_rRNA/run_rnammer.sh bac /absolute/path/to/bacterial_fna /absolute/path/to/new_bac_rnammer /absolute/path/to/rnammer
bash scripts/16s_rRNA/run_rnammer.sh arc /absolute/path/to/archaeal_fna /absolute/path/to/new_arc_rnammer /absolute/path/to/rnammer
python scripts/16s_rRNA/process_rnammer.py --genomes /absolute/path/to/bacterial_fna --gff-dir /absolute/path/to/new_bac_rnammer --out /absolute/path/to/new_bac_16s
python scripts/16s_rRNA/process_rnammer.py --genomes /absolute/path/to/archaeal_fna --gff-dir /absolute/path/to/new_arc_rnammer --out /absolute/path/to/new_arc_16s
```

RNAmmer uses `-S bac` or `-S arc`, `-m ssu`, `-multi`, and outputs per-genome XML, GFF and FASTA. `process_rnammer.py` obtains each sequence from GFF genomic coordinates (reverse-complements minus-strand calls), retains 1450–1700 nt sequences and selects the maximum RNAmmer GFF score per genome. Outputs: `full_length_16s.fna`, `representative_16s.fna`, `16s_calls.tsv` (genomic coordinates, lengths and scores), and `16s_recovery_by_genome.tsv`. Provide `--taxonomy /absolute/path/to/data/final_clean_taxonomy_30495_UPDATED.tsv` to additionally generate `16s_recovery_by_phylum.tsv`. Verify GFF-to-FASTA coordinate consistency, lengths, and expected genome-level counts.

## 3. 16S conservation

The conservation scripts provided by collaborators used MAFFT `--auto --thread 16`, removed alignment columns with more than 50% non-ACGTU characters, calculated Shannon information `max(0,2-H)`, and smoothed with a 12-position moving average. Their code treated `T` and `U` as distinct accepted characters. `conservation.py` preserves that convention by default, with `--normalize-u` providing an optional four-state check.

```bash
mafft --auto --thread 16 /absolute/path/to/full_length_16s.fna > /absolute/path/to/aligned_16s.fna
python scripts/16s_rRNA/conservation.py --alignment /absolute/path/to/aligned_16s.fna --output-prefix /absolute/path/to/bacteria_conservation
```

Outputs: one tabular information profile (`.tsv`) and one final figure (`.pdf`). Original plotted V1–V9 coordinate bands were approximate and should not be interpreted as experimentally mapped positions on an unreferenced alignment.

## 4. VSEARCH v2.31.0

The manuscript describes `--usearch_global` with `--self`, `--maxaccepts 0`, `--maxrejects 0`, `--id 0.0` and `--iddef 2`. The alternative `--allpairs_global --acceptall --qmask none` mode yields a different comparison set. Use the mode corresponding to the intended analysis and check self-hits and pair counts.

```bash
bash scripts/16s_rRNA/run_vsearch.sh manuscript /absolute/path/to/representative_16s.fna /absolute/path/to/vsearch_blast6.tsv
bash scripts/16s_rRNA/run_vsearch.sh allpairs /absolute/path/to/full_length_16s.fna /absolute/path/to/vsearch_allpairs.tsv
```

## 5. Genome-wide scores and statistical analysis

Generate FastANI v1.33, Mash v2.3 (`k=21`, sketch size 1000), and FastAAI v0.1.17 scores following [`30K_FastANI_Mash_FastAAI_Commands.md`](30K_FastANI_Mash_FastAAI_Commands.md). Stage the scores and reference taxonomy for downstream statistical analysis. See `scripts/statistical_analysis/README.md` for required filenames and execution order.

Check expected counts: 30,495 accessions, 16,402 same-species pairs, and 1,754,127 same-genus/different-species pairs.

