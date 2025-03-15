#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Jk_allele_count
#SBATCH --array=1-109
#SBATCH --output=logs/Purging/Jackknife/%x_%j_%a.out
#SBATCH --error=logs/Purging/Jackknife/%x_%j_%a.err
#SBATCH -p shared
#SBATCH -t 3:00:00
#SBATCH --mem=500MB  
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Load required module
module load bcftools/1.20

# Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife

# Get current block coordinates from BED file
current_block=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${savedir}/renamed_blocks.bed | awk '{print $1 ":" $2+1 "-" $3}')

# Process VCF while excluding current block
bcftools view -t ^"${current_block}" "${datadir}/All_Origins_annotated.vcf.gz" | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN[\t%SAMPLE\t%GT]\n' | \
awk -v block=${SLURM_ARRAY_TASK_ID} '
BEGIN {
    OFS = "\t";
    print "Sample", "HIGH", "MODERATE", "LOW", "MODIFIER", "Total_Alleles";
}
{
    chrom = $1
    pos = $2
    ref = $3
    alt = $4
    
    split($5, ann_list, ",");
    for (a in ann_list) {
        split(ann_list[a], fields, "|");
        impact = fields[3];
        effect_allele = substr(fields[1], 1, 1);

        if (effect_allele == ref) continue;

        for (i = 6; i <= NF; i += 2) {
            sample = $i
            gt = $(i+1)

            if (gt ~ /^\./) continue

            split(gt, alleles, /[\/|]/)
            total = 0
            alt_count = 0
            for (a in alleles) {
                if (alleles[a] != ".") total++
                if (alleles[a] == "1") alt_count++
            }
            total_alleles[sample] += total

            if (impact == "HIGH") counts[sample]["HIGH"] += alt_count
            else if (impact == "MODERATE") counts[sample]["MODERATE"] += alt_count
            else if (impact == "LOW") counts[sample]["LOW"] += alt_count
            else if (impact == "MODIFIER") counts[sample]["MODIFIER"] += alt_count
        }
    }
}
END {
    for (sample in counts) {
        printf "%s\t%d\t%d\t%d\t%d\t%d\n",
            sample,
            counts[sample]["HIGH"] + 0,
            counts[sample]["MODERATE"] + 0,
            counts[sample]["LOW"] + 0,
            counts[sample]["MODIFIER"] + 0,
            total_alleles[sample] + 0;
    }
}' > "${savedir}/Jakknife_block${SLURM_ARRAY_TASK_ID}.tsv"

echo "Results saved to ${savedir}/Jakknife_block${SLURM_ARRAY_TASK_ID}.tsv"