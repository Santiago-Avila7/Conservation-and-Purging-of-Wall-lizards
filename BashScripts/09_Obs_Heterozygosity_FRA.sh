#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J Obs_Het 
#SBATCH --output=logs/Gene_diversity/%x_%j.out
#SBATCH --error=logs/Gene_diversity/%x_%j.err
#SBATCH -p shared
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Count heterozygous sites per sample in French-origin samples 

# Load bcftools module
module load bcftools/1.20

# Define directories
datadir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data"
mkdir -p "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity"
savedir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity"

# Define file paths
VCF_GZ="${datadir}/All_French_final.vcf.gz"
OUTPUT="${savedir}/All_French_het.tsv"

# Process each individual listed in the individuals.txt file
echo "Processing VCF..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' "$VCF_GZ" | awk '
BEGIN { 
    OFS="\t"; 
    print "Sample", "Heterozygous_Count", "Total_Variant_Sites", "Missing_Genotypes";
}
{
    has_alt = ($4 != ".") ? 1 : 0;  # Ensure the site has at least one ALT allele

    if (has_alt) {
        for (i = 5; i <= NF; i++) {  
            split($i, sample_gt, "=");  # Extract sample name and genotype
            sample = sample_gt[1];
            gt = sample_gt[2];

            # Initialize sample counters if not already done
            if (!(sample in het_count)) {
                het_count[sample] = 0;
                total_variants[sample] = 0;
                missing_count[sample] = 0;
            }

            total_variants[sample]++;  # Count this site as a variant site for the sample

            if (gt == "./." || gt == ".|.") {
                missing_count[sample]++;
            } else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
                het_count[sample]++;
            }
        }
    }
}
END {
    for (sample in het_count) {
        print sample, het_count[sample], total_variants[sample], missing_count[sample];
    }
}' > "$OUTPUT"

echo "Results saved to $OUTPUT"
