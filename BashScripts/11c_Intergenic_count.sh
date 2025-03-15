#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Intergenic_count 
#SBATCH --output=logs/Purging/%x_%j.out
#SBATCH --error=logs/Purging/%x_%j.err
#SBATCH -p shared
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
##SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se


# Load required module
module load bcftools/1.20

# Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts

# Filter importat positions and pipe for allele count. 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE\t%GT]\n' ${datadir}/All_intergenic.vcf.gz | \
awk '
BEGIN {
    FS="\t";  # Set field separator to tab
    OFS="\t"; # Set output field separator to tab
    print "Sample", "Total_Alternates", "Total_Variants";  # Print header
}
{
    # Loop through each sample-genotype pair
    for (i = 5; i <= NF; i += 2) {
        sample = $i;      # Extract sample name
        gt = $(i + 1);    # Extract genotype

        # Skip missing genotypes (./. or .|.)
        if (gt ~ /^\./) continue;

        # Count total variants (non-missing genotypes)
        total_variants[sample]++;

        # Count alternate alleles
        split(gt, alleles, /[\/|]/);  # Split genotype into alleles
        for (a in alleles) {
            if (alleles[a] == "1") {  # Count alternate alleles (1)
                total_alternates[sample]++;
            }
        }
    }
}
END {
    # Print results for each sample
    for (sample in total_variants) {
        print sample, total_alternates[sample] + 0, total_variants[sample] + 0;
    }
}' > ${datadir}/intergenic_allele_counts.tsv