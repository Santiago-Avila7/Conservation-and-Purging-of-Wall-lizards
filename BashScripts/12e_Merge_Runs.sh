#!/bin/bash

# Output directory and file
dir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/French_Italian_Native"
output_file="${dir}/All_Natives_jk.tsv"

# Write header (matches the new script's output format)
echo -e "Jackknife\tOrigin\tAlternate_HIGH\tTotal_HIGH\tAlternate_MODERATE\tTotal_MODERATE\tAlternate_LOW\tTotal_LOW\tAlternate_MODIFIER\tTotal_MODIFIER\tAlternate_INTERGENIC\tTotal_INTERGENIC" > "$output_file"

# Process each jackknife file
for f in "${dir}"/*_results.tsv; do
    # Extract jackknife ID from filename (e.g., jackknife_42_results.tsv → 42)
    jk=$(basename "$f" | grep -oP '(?<=jackknife_)[0-9]+(?=_results)')
    
    # Append data while:
    # 1. Skipping the header (NR > 1)
    # 2. Adding the jackknife ID as first column
    awk -v jk="$jk" 'NR > 1 {print jk "\t" $0}' "$f" >> "$output_file"
done

# Verify output
echo "Merged results saved to: $output_file"
wc -l "$output_file"