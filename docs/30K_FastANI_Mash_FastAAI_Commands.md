# FastANI, Mash, and FastAAI: commands and parallel execution

**Project:** 30,495 bacterial and archaeal type-material assemblies  
**Purpose:** Generate FastANI, Mash, and FastAAI scores using genome lists, individual pairs, or parallel batches.

## 1. Software, inputs, and execution order

| Software | Manuscript version | Input | Core output | Parallel/list option |
|---|---|---|---|---|
| FastANI | 1.33 | Assembly FASTA files | Query, reference, ANI, matching fragments, query fragments | `--ql`, `--rl`, `-t` |
| Mash | 2.3 | Assembly FASTA files or sketches | Reference, query, distance, p-value, shared hashes | `mash sketch -l`, `-p`; `mash dist -p` |
| FastAAI | 0.1.17 | One nucleotide FASTA file per genome | Raw mean tetramer Jaccard, shared SCPs, estimated AAI | `build_db --genomes` accepts a **path-list file**; `--threads` |

Run order: **freeze accession/FASTA manifest → choose intended comparison pairs → generate scores → verify pair coverage and software versions → select manuscript SS and SGDS comparisons → downstream statistical scripts.** The frozen input accession list is `data/gca_accessions_30495.txt` in the project repository. Each program's outputs should retain identifiers that can be mapped back to the exact GCA accession.version.

**A full 30,495 × 30,495 comparison is not required for the species-discrimination dataset.** There are 464,957,265 distinct unordered pairs in that exhaustive universe, whereas the reported species discrimination uses specifically defined same-species (SS) and same-genus/different-species (SGDS) pairs. A list-based *Cartesian product* is not equivalent to an arbitrary selected-pair list.

### Prepare a real FASTA manifest (edit the absolute paths)

```bash
#!/usr/bin/env bash
set -euo pipefail
GENOMES=/absolute/path/to/30k/genomes
WORK=/absolute/path/to/30k/work
mkdir -p "$WORK/lists" "$WORK/results" "$WORK/logs"

# One assembly FASTA file per genome, with stable GCA accession.version in filename.
find "$GENOMES" -type f \( -name '*.fna' -o -name '*.fna.gz' \) -print \
  | LC_ALL=C sort > "$WORK/lists/genome_paths.txt"
wc -l "$WORK/lists/genome_paths.txt"

# A ten-real-genome smoke-test subset, not the study's final pair selection.
head -n 10 "$WORK/lists/genome_paths.txt" > "$WORK/lists/pilot10.txt"

# Capture local versions and tool help before production use.
command -v fastANI
command -v mash
command -v fastaai
mash --version
python -m pip show FastAAI
fastANI -h > "$WORK/logs/fastani_help.txt" 2>&1 || true
fastaai build_db -h > "$WORK/logs/fastaai_build_help.txt" 2>&1 || true
fastaai db_query -h > "$WORK/logs/fastaai_query_help.txt" 2>&1 || true
```

Use the frozen accession-to-FASTA mapping. The pilot list is for testing commands. Check installed software versions before execution.

## 2. FastANI v1.33

### Single genome pair

```bash
fastANI \
  -q /absolute/path/to/30k/genomes/GCA_QUERY.fna \
  -r /absolute/path/to/30k/genomes/GCA_REFERENCE.fna \
  -o /absolute/path/to/30k/work/results/fastani_one_pair.tsv
```

### Query and reference genome lists

Each list must contain **one absolute path to an assembly FASTA per line**. `--ql` and `--rl` calculate the query × reference Cartesian product. `-t` is internal multithreading.

```bash
fastANI \
  --ql /absolute/path/to/30k/work/lists/query_genomes.txt \
  --rl /absolute/path/to/30k/work/lists/reference_genomes.txt \
  -t 16 \
  -o /absolute/path/to/30k/work/results/fastani_batch.tsv
```

**List-based command:** `fastANI --ql genomes.txt --rl genomes.txt -t 32 -o fastani_all_vs_all.tsv`. FastANI v1.33 defaults include a 3,000 bp fragment length and minimum matching fraction of 0.2.

FastANI's tab-separated output has five fields:

```text
query_fasta    reference_fasta    ANI_percent    matched_fragments    total_query_fragments
```

Compute **query-direction alignment fraction** as `matched_fragments / total_query_fragments`. Preserve all five source columns for **both SS and SGDS**. Missing FastANI output for an intended pair must be tracked as non-recovery, never filled with zero or ANI=0. FastANI can be direction-dependent, so record how A→B, B→A and duplicate records were handled.

### Optional bounded GNU Parallel execution over query batches

First create nonoverlapping query lists, each containing absolute FASTA paths, in `/absolute/path/to/30k/work/lists/query_batches/`. Use the **same frozen reference list** for all batches. This evaluates the defined query batches against those references. It is **not** a runner for arbitrary selected pairs.

```bash
#!/usr/bin/env bash
set -euo pipefail
WORK=/absolute/path/to/30k/work
mkdir -p "$WORK/results/fastani_batches"

find "$WORK/lists/query_batches" -maxdepth 1 -type f -name 'queries_*.txt' -print0 \
  | sort -z \
  | parallel -0 --jobs 4 --halt soon,fail=1 \
      'fastANI --ql {} --rl /absolute/path/to/30k/work/lists/reference_genomes.txt -t 4 -o /absolute/path/to/30k/work/results/fastani_batches/{/.}.tsv'
```

Four simultaneous jobs × four FastANI threads each can require approximately 16 CPU cores, plus memory for each process. Do not launch millions of individual single-pair processes as a default production strategy. Validate the number of expected query–reference records, including non-recovered pairs, before merging results.

## 3. Mash v2.3

### Sketch assemblies from a genome list

```bash
mash sketch \
  -k 21 -s 1000 -p 16 \
  -l /absolute/path/to/30k/work/lists/pilot10.txt \
  -o /absolute/path/to/30k/work/results/pilot10_mash
```

This creates `pilot10_mash.msh` with one **whole-genome** sketch per FASTA file. For multi-contig assemblies, **do not add `-i`**: `-i` sketches each sequence/contig separately, which changes the unit of analysis. The manuscript reports `k=21`, sketch size `1000`.

### Calculate pairwise sketch distances

```bash
mash dist -p 16 \
  /absolute/path/to/30k/work/results/pilot10_mash.msh \
  /absolute/path/to/30k/work/results/pilot10_mash.msh \
  > /absolute/path/to/30k/work/results/pilot10_mash_dist.tsv
```

Mash emits `reference_ID`, `query_ID`, `distance`, `p_value`, `shared_hashes`. A same-sketch × same-sketch run includes self-comparisons and both pair directions. Subset to the exact intended comparison pairs and deduplicate by a documented rule before manuscript analysis. **Mash distance=1 with no detectable shared sketch hashes is a valid numeric result, not an absent/failed output**. Preserve it distinctly so downstream sensitivity choices remain explicit.

`mash sketch -l` supports genome-path lists and `-p` sets threads. For `mash dist`, `-p` sets threads. `mash dist -l` instead treats **query arguments** as files containing query-path lists, not as the same type of input as `mash sketch -l`. The commands above use precomputed sketches to avoid that ambiguity.

## 4. FastAAI v0.1.17

### Build a database from genome paths

FastAAI `--genomes` accepts a directory, **a text file with one genome FASTA path per line**, or comma-separated file paths. The text-file option is preferable for controlling the accession set. Each genome must occupy one nucleotide FASTA file.

```bash
fastaai build_db \
  --genomes /absolute/path/to/30k/work/lists/pilot10.txt \
  --threads 8 \
  --output /absolute/path/to/30k/work/results/fastaai_pilot \
  --database pilot10.db \
  --compress
```

The database should be created at:

```text
/absolute/path/to/30k/work/results/fastaai_pilot/database/pilot10.db
```

### Query a database against itself

```bash
fastaai db_query \
  --query /absolute/path/to/30k/work/results/fastaai_pilot/database/pilot10.db \
  --target /absolute/path/to/30k/work/results/fastaai_pilot/database/pilot10.db \
  --threads 8 \
  --output_style tsv \
  --output /absolute/path/to/30k/work/results/fastaai_pilot
```

FastAAI's tabular output contains `query`, `target`, **`avg_jacc_sim` (raw mean tetramer Jaccard)**, `jacc_SD`, `num_shared_SCPs`, `poss_shared_SCPs`, and `AAI_estimate`. The paper uses **raw Jaccard**, not the transformed estimated AAI column. Retain the shared-SCP counts and missing/no-shared-SCP cases. Verify the installed FastAAI v0.1.17 CLI with `fastaai build_db -h` and `fastaai db_query -h`.

FastAAI `--threads` parallelizes processing within the program. For separate taxon-specific groups, independent database builds/queries may be scheduled as separate jobs when comparison-pair selection and database membership are preserved. Do not assume that arbitrary disjoint groups are equivalent to cross-group or domain-wide comparisons.

## 5. Select intended pairs and validate output

Prepare the intended-pair table, **not merely a complete list of genomes**. Suggested TSV columns:

```text
gca_query  gca_reference  taxonomic_relation  method_expected
```

Recommended reconciliation checks:

1. Confirm the 30,495 accession.version set and fixed SS/SGDS pair lists; do not refresh taxonomy during the analysis.
2. Map all FASTA paths and FastAAI/Mash identifiers to unambiguous GCA accession.version values.
3. For each method, report expected pair count, observed numeric result count, non-recovered count, duplicated/reversed pairs, and records with invalid formats.
4. For FastANI, retain directional ANI and both fragment counts, and calculate alignment fraction separately for SS and SGDS.
5. For Mash, retain distance=1 separately from absent results. Do not confuse sketch non-detection with program failure.
6. For FastAAI, use **raw `avg_jacc_sim`**, not `AAI_estimate`.
7. Record program versions, input checksums, commands, thread counts, output checksums, and filtering or score aggregation before running `scripts/statistical_analysis/final_species_discrimination/scripts/prepare_inputs.py`.

## 6. Software documentation

- FastANI: <https://github.com/ParBLiSS/FastANI> — list options `--ql`, `--rl`, `-t`, five output fields, and query/reference asymmetry.
- Mash: <https://github.com/marbl/Mash> and <https://mash.readthedocs.io/en/latest/> — `sketch -l`, `-k`, `-s`, `-p`, `dist`, whole-file versus individual-contig sketches.
- FastAAI: <https://github.com/cruizperez/FastAAI> — `build_db`, `db_query`, list-of-genomes input, threaded execution, and raw Jaccard output.
- Version record for FastAAI v0.1.17: <https://pypi.org/project/FastAAI/0.1.17/>.
