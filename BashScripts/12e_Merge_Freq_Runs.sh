#!/bin/bash

INPUT_DIR="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/All_Italian"
OUTPUT_FILE="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/All_Italian/Italian_Freq_results.tsv"

# Write header with Jackknife_ID column
echo -n "Jackknife_ID" > "$OUTPUT_FILE"
head -n 1 $(ls ${INPUT_DIR}/jackknife_*_relative.tsv | head -1) >> "$OUTPUT_FILE"

# Append each file's data with ID
for file in ${INPUT_DIR}/jackknife_*_relative.tsv; do
    jk_id=$(basename "$file" | cut -d'_' -f2)
    tail -n +2 "$file" | awk -v id="$jk_id" '{print id "\t" $0}' >> "$OUTPUT_FILE"
done