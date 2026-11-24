#!/usr/bin/env bash
# Usage: bash fetch_ncbi_taxids.sh /absolute/manifest.txt /absolute/new_output_dir
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'Usage: bash fetch_ncbi_taxids.sh MANIFEST OUTPUT_DIR' >&2; exit 2; }
manifest=$1; out=$2
[[ -s "$manifest" ]] || { echo 'Missing accession manifest' >&2; exit 1; }
[[ ! -e "$out" ]] || { echo 'Output directory exists; refusing to overwrite' >&2; exit 1; }
command -v datasets >/dev/null || { echo 'NCBI datasets CLI not in PATH' >&2; exit 1; }
command -v jq >/dev/null || { echo 'jq not in PATH' >&2; exit 1; }
mkdir -p "$out"
datasets summary genome accession --inputfile "$manifest" --as-json-lines > "$out/assembly_summary.jsonl"
jq -r '[.accession, .organism.tax_id] | @tsv' "$out/assembly_summary.jsonl" > "$out/assembly_taxids.tsv"
cut -f2 "$out/assembly_taxids.tsv" | grep -v '^null$' | sort -u > "$out/taxids.txt"
[[ -s "$out/taxids.txt" ]] || { echo 'TaxID extraction returned no records' >&2; exit 1; }
datasets summary taxonomy taxon --inputfile "$out/taxids.txt" --as-json-lines > "$out/taxonomy.jsonl"
printf 'Assemblies: %s; non-null distinct TaxIDs: %s\n' "$(wc -l < "$out/assembly_taxids.tsv")" "$(wc -l < "$out/taxids.txt")"
