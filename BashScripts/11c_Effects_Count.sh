#!/bin/bash

# Count the number of sites based on the impact category

mamba activate snpeff #Activate before running 

zcat /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/All_Origins_annotated.vcf.gz | \
SnpSift extractFields -e "." - CHROM POS ANN[0].IMPACT | \ 
sort -u | \
cut -f 3 | \
sort | \
uniq -c | \
awk '{print $2 "\t" $1}' > /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Effect_Counts.tsv

# The first line in the output is a register of "ANN[0].IMPACT" that correspond to the header. It is not a real annotation in a variant