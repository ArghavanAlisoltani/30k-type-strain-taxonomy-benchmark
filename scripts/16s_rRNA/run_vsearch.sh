#!/usr/bin/env bash
# Usage: bash run_vsearch.sh manuscript|allpairs INPUT_FASTA OUTPUT_TSV
set -euo pipefail
if [[ $# -ne 3 || ( "$1" != "manuscript" && "$1" != "allpairs" ) ]]; then
  echo 'Usage: bash run_vsearch.sh manuscript|allpairs INPUT_FASTA OUTPUT_TSV' >&2
  exit 2
fi
mode=$1; input=$2; output=$3
[[ -s "$input" ]] || { echo 'Missing input FASTA' >&2; exit 1; }
[[ ! -e "$output" ]] || { echo 'Refusing to overwrite output' >&2; exit 1; }
command -v vsearch >/dev/null || { echo 'vsearch not found in PATH' >&2; exit 1; }
if [[ "$mode" == manuscript ]]; then
  # The manuscript reports --usearch_global and --self; confirm historical command/log.
  vsearch --usearch_global "$input" --db "$input" --maxaccepts 0 --maxrejects 0 \
    --self --id 0.0 --iddef 2 --blast6out "$output"
else
  # The newly supplied methods note uses unique unordered all-versus-all pairs.
  vsearch --allpairs_global "$input" --acceptall --qmask none --iddef 2 \
    --userout "$output" --userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi
fi
