# Genome comparison and taxonomy utilities

## Software

- Python 3.9+ (standard library only)
- FastANI 1.33 for the blockwise runner
- Mash 2.3 for triangular distance matrices

## 1. Check the taxonomy dictionary

Use an input `GenID={...}` dictionary with the frozen accession manifest. The 30,495-genome dataset excludes `GCA_041504555.1` and `GCA_041505935.1`.

```bash
python /absolute/path/to/repo/scripts/taxonomy/audit_dictionary.py \
  --dictionary /absolute/path/to/dictionary.py \
  --manifest /absolute/path/to/repo/data/gca_accessions_30495.txt \
  --reference-tsv /absolute/path/to/repo/data/final_clean_taxonomy_30495_UPDATED.tsv
```

The script reports accession and rank-label differences without writing files. Add `--export-tsv /absolute/path/to/new_taxonomy_subset.tsv` to write a new TSV limited to the frozen accession list. The output is a new file; the input reference TSV is unchanged.

## 2. FastANI in blocks

`--genomes` accepts an absolute directory of `.fna`, `.fa`, `.fasta` (optionally `.gz`) files or an absolute text file containing one absolute genome FASTA path per line.

```bash
# Plan only: show blocks, pairs, and CPU request without computing anything.
python /absolute/path/to/repo/scripts/fastani/fastani_blockwise.py \
  --genomes /absolute/path/to/query_genomes.txt \
  --out /absolute/path/to/new_fastani_run \
  --batch-size 500 --threads 16 --jobs 2

# Run FastANI on all unordered pairs in one genome set.
python /absolute/path/to/repo/scripts/fastani/fastani_blockwise.py \
  --genomes /absolute/path/to/query_genomes.txt \
  --out /absolute/path/to/new_fastani_run \
  --batch-size 500 --threads 16 --jobs 2 --execute
```

Use `--references /absolute/path/to/reference_genomes.txt` for a **query × reference** comparison instead. `--jobs` controls concurrent FastANI processes; `--threads` controls threads per process. Reserve at least `jobs × threads` CPU cores and enough memory.

Defaults: `--frag-len 3000`, `--min-fraction 0.2`. The runner uses FastANI's `--ql`/`--rl` options internally. Completed blocks have `.done` markers; failed jobs are never marked complete. Output files:

- `run_config.json` — input-list hashes and parameters
- `lists/*.txt`, `blocks/*.tsv`, `blocks/*.done`, `logs/*.log` — resumable jobs
- `fastani_pairs.tsv` — **no header**; query, reference, ANI %, matched fragments, total query fragments
- `fastani_summary.json` — possible pairs, numeric pairs, pairs without numeric output

Genome path identifiers are preserved. For diagonal blocks, the script retains the lexicographically ordered query–reference direction when both directions are available. Fragment counts are **directional** and are never associated with reversed identifiers. Query-direction alignment fraction is `matched_fragments / total_query_fragments`. An absent FastANI result is not ANI=0.

For a query set of independent non-type genomes, use separate, biologically relevant type/reference groups. The list mode is a Cartesian product, not an exact arbitrary-pair selector.

To merge completed blocks without rerunning FastANI:

```bash
python /absolute/path/to/repo/scripts/fastani/fastani_blockwise.py \
  --genomes /absolute/path/to/query_genomes.txt \
  --out /absolute/path/to/new_fastani_run \
  --batch-size 500 --threads 16 --jobs 2 --merge-only
```

## 3. Mash triangular matrix to pairs

Use only for **`mash triangle`** lower-triangular relaxed-PHYLIP output. The first line declares the number of genomes, followed by rows containing a genome ID and its distances to earlier genomes. Genome IDs must contain no whitespace.

```bash
python /absolute/path/to/repo/scripts/mash/mash_triangle_to_pairs.py \
  /absolute/path/to/mash_triangle.txt \
  /absolute/path/to/mash_pairs.tsv
```

Output is **no header**, three tab-separated columns: `genome_a`, `genome_b`, `Mash_distance`. The converter validates count, row widths, IDs, and distances before creating the output. The standard five-column `mash dist` output is already an edge list and does **not** require conversion.

For the Mash sketching/distance and FastAAI building/query commands, see `docs/30K_FastANI_Mash_FastAAI_Commands.md`.

## 4. Tests

```bash
python -m unittest discover -s /absolute/path/to/repo/tests -p 'test_score_generation.py' -v
```

Unit tests cover parsing, directionality, incomplete jobs, and accession filtering.
