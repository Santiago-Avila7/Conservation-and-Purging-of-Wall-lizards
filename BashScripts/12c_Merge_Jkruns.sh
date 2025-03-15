#!/bin/bash

# Merge tvs per Jakknie
mergedir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jakknife"

# Create header with new column
echo -e "Jakknife\tSample\tHIGH\tMODERATE\tLOW\tMODIFIER\tTotal_Alleles" > "${mergedir}/merged_jackknife.tsv"

# Process files in numerical order
for block in {1..109}; do
    file="${mergedir}/Jakknife_block${block}.tsv"
    if [[ -f "$file" ]]; then
        # Add block number column and skip header
        awk -v blk="$block" 'NR>1 {print blk "\t" $0}' "$file" >> "${mergedir}/merged_jackknife.tsv"
    fi
done