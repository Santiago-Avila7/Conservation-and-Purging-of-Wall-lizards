#!/bin/bash

INPUT_DIR="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/IntroducedVsNative/Italian"
OUTPUT_FILE="${INPUT_DIR}/All_Italian_freq.tsv"

# Get header from first file and add Jackknife_ID column with TAB separator
FIRST_FILE=$(ls ${INPUT_DIR}/jackknife_*_relative.tsv | head -1)
HEADER=$(head -n 1 "$FIRST_FILE")
echo -e "Jackknife_ID\t${HEADER}" > "$OUTPUT_FILE"

# Process each file
for file in ${INPUT_DIR}/jackknife_*_relative.tsv; do
    # Extract jackknife ID from filename (e.g., extract "5" from "jackknife_5_relative.tsv")
    jk_id=$(basename "$file" | cut -d'_' -f2)
    
    # Get data line (line 2) and prepend ID with TAB separator
    DATA_LINE=$(sed -n '2p' "$file")
    echo -e "${jk_id}\t${DATA_LINE}" >> "$OUTPUT_FILE"
done

echo "Merged $(ls ${INPUT_DIR}/jackknife_*_relative.tsv | wc -l) files into $OUTPUT_FILE"