#!/bin/bash

# Output file
dir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife2

# Write header once
echo -e "Jackknife\tOrigin\tHIGH\tMODERATE\tLOW\tMODIFIER\tTotal_Alleles\tIntergenic_Alternate\tIntergenic_Total" > "${dir}//all_results.tsv"

# Process each jackknife file
for f in ${dir}/*.tsv; do
    # Extract jackknife ID from filename
    jk=$(basename "$f" | grep -oP '[0-9]+(?=_)')
    
    # Append data to output file (skip header in input files)
    awk -v jk="$jk" 'NR > 1 {print jk "\t" $0}' "$f" >> "${dir}//all_results.tsv"
done