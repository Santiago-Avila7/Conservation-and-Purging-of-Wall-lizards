#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Countcorrection_Pmuralis
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.err
#SBATCH -p shared
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Corrects the outgroup nucleotides counts since est-sfs allows just one nucleotide per outgroup

# Files
INPUT="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pmuralis_counts.txt"
OUTPUT="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pmuralis_counts_corr.txt"


awk '
function is_homozygous(a,b,c,d) {
    return ((a==2 && b==0 && c==0 && d==0) || 
            (a==0 && b==2 && c==0 && d==0) || 
            (a==0 && b==0 && c==2 && d==0) || 
            (a==0 && b==0 && c==0 && d==2))
}
function is_heterozygous(a,b,c,d) {
    return ((a==1 && b==1 && c==0 && d==0) || 
            (a==1 && b==0 && c==1 && d==0) || 
            (a==1 && b==0 && c==0 && d==1) || 
            (a==0 && b==1 && c==1 && d==0) || 
            (a==0 && b==1 && c==0 && d==1) || 
            (a==0 && b==0 && c==1 && d==1))
}
{
    printf "%s", $1
    for (i=2; i<=4; i++) {
        split($i, counts, ",")
        a=counts[1]+0; b=counts[2]+0; c=counts[3]+0; d=counts[4]+0
        if (is_heterozygous(a,b,c,d)) het_outgroups++
        og[i]=sprintf("%d,%d,%d,%d",a,b,c,d)
    }
    for (i=2; i<=4; i++) {
        split(og[i], counts, ",")
        a=counts[1]+0; b=counts[2]+0; c=counts[3]+0; d=counts[4]+0
        if (is_homozygous(a,b,c,d)) {
            if (a==2) printf "\t1,0,0,0"
            else if (b==2) printf "\t0,1,0,0"
            else if (c==2) printf "\t0,0,1,0"
            else if (d==2) printf "\t0,0,0,1"
        } else if (is_heterozygous(a,b,c,d)) {
            if (het_outgroups > 1) printf "\t0,0,0,0"
            else printf "\t%s", og[i]
        } else {
            printf "\t0,0,0,0"
        }
    }
    het_outgroups=0
    printf "\n"
}' "$INPUT" > "$OUTPUT"

echo "Processed counts written to: $OUTPUT"
