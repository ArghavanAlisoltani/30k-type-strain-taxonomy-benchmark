#!/bin/bash

# Ensure RNAmmer is in your PATH or update the variable below
RNAMMER_PATH="/home/ranabir/Prof_Dave/rnammer-1.2.Unix/rnammer"
mkdir -p results

for file in *.fna; do
    base="${file%.fna}"
    
    echo "Processing $file..."
           
    # Run RNAmmer and output to FASTA (.fna) (-f) lsu,ssu,tsu
    perl "$RNAMMER_PATH" -S bac -m ssu  \
        -xml "results/${base}.xml" \
        -gff "results/${base}.gff" \
        -f "results/${base}_16s_rRNAs.fasta" \
        -multi \
        < "$file"

    #echo "Extracting 18s (SSU) for $base..."
    # Extracts the 18s sequences into individual files
    #awk '/^>/ {print_seq=(index($0,"18s_rRNA") > 0)} print_seq' "results/${base}_all_rRNAs.fasta" > "results/${base}_18s.fasta"
done

# --- NEW: Concatenation Step ---
echo "Merging all 16s sequences into master file..."
cat results/*_16s_rRNAs.fasta > results/master_16s_rRNA.fasta

echo "Done! Final master file: master_16s_rRNA.fna"
