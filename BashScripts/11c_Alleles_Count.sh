#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Ind_Effects 
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
# Define directories and genome name
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts

# Define file paths
VCF_GZ="${datadir}/All_Origins_annotated.vcf.gz"
OUTPUT="${savedir}/alleles_per_sample.tsv"

# Extract annotations & genotypes while **keeping correct sample names**
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN[\t%SAMPLE\t%GT]\n' "$VCF_GZ" | awk '
BEGIN {
    OFS = "\t";
    print "Sample", "HIGH", "MODERATE", "LOW", "MODIFIER", "Total_Alleles", "Derived_Alleles", "Effect_Allele_is_REF";
}
{
    split($5, ann_list, ",");  # Split multi-annotations by commas
    for (a in ann_list) {
        split(ann_list[a], fields, "|");  # Split SnpEff annotation by "|"
        impact = fields[3];  # Impact type
        effect_allele = substr(fields[1], 1, 1);  # First nucleotide of ANN field

        # Ensure the effect allele matches ALT (not REF)
        if (effect_allele != $4) continue;

        for (i = 6; i <= NF; i += 2) {  # Iterate through sample-genotype pairs
            sample = $i;    # Extract sample name
            gt = $(i + 1);  # Extract genotype

            if (gt != "./." && gt != ".|.") {  # Ignore missing data
                total_alleles[sample] += 2;  # Each individual has 2 alleles per site

                # Count derived alleles (ALT)
                if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
                    derived_alleles[sample] += 1;  # Heterozygous contributes 1 derived allele
                }
                if (gt == "1/1" || gt == "1|1") {
                    derived_alleles[sample] += 2;  # Homozygous ALT contributes 2 derived alleles
                }

                # Count cases where the effect allele is actually REF
                if (effect_allele == $3) {
                    effect_allele_is_ref[sample]++;
                }

                # Count by impact category
                if (impact == "HIGH") counts[sample]["HIGH"]++;
                else if (impact == "MODERATE") counts[sample]["MODERATE"]++;
                else if (impact == "LOW") counts[sample]["LOW"]++;
                else if (impact == "MODIFIER") counts[sample]["MODIFIER"]++;
            }
        }
    }
}
END {
    for (sample in counts) {
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            sample,
            counts[sample]["HIGH"] + 0,
            counts[sample]["MODERATE"] + 0,
            counts[sample]["LOW"] + 0,
            counts[sample]["MODIFIER"] + 0,
            total_alleles[sample] + 0,
            derived_alleles[sample] + 0,
            effect_allele_is_ref[sample] + 0;
    }
}' > "$OUTPUT"

echo "Results saved to $OUTPUT"
