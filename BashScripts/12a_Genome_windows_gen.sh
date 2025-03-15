#!/bin/bash

# Directories 
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jakknife
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jakknife 

# Module 
module load bedtools/2.31.0

# Filter autosomes (adjust regex to match your naming convention)
grep -E '^CM[0-9]{6}\.[0-9]+' ${datadir}/GCA_004329235.1_PodMur_1.0_genomic.fna.fai > ${savedir}/autosomes.fai

# Create genome file for bedtools
awk '{print $1 "\t" $2}' ${savedir}/autosomes.fai > ${savedir}/autosomes.genome
sed -i '/CM014761\.1/d' ${savedir}/autosomes.genome

echo "Genome file generated" 

# Calculate total autosome length
total_length=$(awk '{sum += $2} END {print sum}' ${savedir}/autosomes.genome)

# Split into 100 blocks (each ~1% of total length)
bedtools makewindows \
  -g ${savedir}/autosomes.genome \
  -w $((total_length / 100)) \
  > ${savedir}/autosomes_100_blocks.bed

echo "Genome blocks created"

# Add block IDs
awk '{print $0 "\tblock" NR}' ${savedir}/autosomes_100_blocks.bed > ${savedir}/autosomes_100_blocks_id.bed

echo "Blocks named"

# Rename chromosome per block
sed -f <(awk '{print "s/"$1"/"$2"/"}' Purging_Chr_names.txt) ${savedir}/autosomes_100_blocks_id.bed > ${savedir}/renamed_blocks.bed

echo "Blocks named"