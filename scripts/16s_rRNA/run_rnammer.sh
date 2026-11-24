#!/usr/bin/env bash
# Usage: bash run_rnammer.sh bac|arc /absolute/genomes_dir /absolute/output_dir /absolute/rnammer
set -euo pipefail
if [[ $# -ne 4 || ( "$1" != "bac" && "$1" != "arc" ) ]]; then
  echo 'Usage: bash run_rnammer.sh bac|arc GENOMES_DIR OUTPUT_DIR RNAMMER_PATH' >&2
  exit 2
fi
model=$1; genomes_dir=$2; output_dir=$3; rnammer=$4
[[ -d "$genomes_dir" ]] || { echo 'Genome directory not found' >&2; exit 1; }
[[ -f "$rnammer" ]] || { echo 'RNAmmer executable not found' >&2; exit 1; }
[[ ! -e "$output_dir" ]] || { echo 'Output directory exists; refusing to overwrite' >&2; exit 1; }
mkdir -p "$output_dir"
shopt -s nullglob
files=("$genomes_dir"/*.fna)
(( ${#files[@]} > 0 )) || { echo 'No genome .fna files' >&2; exit 1; }
for file in "${files[@]}"; do
  base=$(basename "$file" .fna)
  perl "$rnammer" -S "$model" -m ssu -multi \
    -xml "$output_dir/$base.xml" -gff "$output_dir/$base.gff" \
    -f "$output_dir/${base}_16s_rRNAs.fasta" < "$file"
done
# Separate per-assembly FASTA/GFF records are authoritative for subsequent extraction.
