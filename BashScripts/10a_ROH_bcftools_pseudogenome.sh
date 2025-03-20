#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J ROH_BCFTools_PSG
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/ROHs/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/ROHs/%x_%j.err
#SBATCH -p shared
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIl
#SBATCH --mail-user=sa5674av-s@student.lu.se


# Load required modules
module load bcftools/1.20

# Define input/output directories
datadir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data"
savedir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity/ROHs"
mkdir -p "$savedir"

# Define input VCFs
french_vcf="${datadir}/All_French_final.vcf.gz"
italian_vcf="${datadir}/All_Italian_final.vcf.gz"

# Define output VCFs with pseudogenome
french_pseudo_vcf="${savedir}/All_French_pseudohom.vcf.gz"
italian_pseudo_vcf="${savedir}/All_Italian_pseudohom.vcf.gz"

# Step 1: Extract VCF header and manually add the sample name
echo "Preparing pseudohom VCFs..."
zcat "$french_vcf" | head -n 1000 | grep "^#" > "${savedir}/header_french.tmp"
zcat "$italian_vcf" | head -n 1000 | grep "^#" > "${savedir}/header_italian.tmp"

# **Manually edit the last line in header_french.tmp and header_italian.tmp to add "pseudohom" as a sample**

# Step 2: Append a fully homozygous (0/0) genotype to all positions
zcat "$french_vcf" | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat "${savedir}/header_french.tmp" - | bgzip > "$french_pseudo_vcf"
zcat "$italian_vcf" | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat "${savedir}/header_italian.tmp" - | bgzip > "$italian_pseudo_vcf"

# Step 3: Index the new pseudo-genome VCFs
bcftools index "$french_pseudo_vcf"
bcftools index "$italian_pseudo_vcf"

# Step 4: Run bcftools roh on modified VCFs
echo "Running ROH analysis..."
bcftools roh --GTs-only 30 --AF-tag AF -Or -o "${savedir}/All_French_ROH_Pseudo_BCFtools" "$french_pseudo_vcf"
bcftools roh --GTs-only 30 --AF-tag AF -Or -o "${savedir}/All_Italian_ROH_Pseudo_BCFtools" "$italian_pseudo_vcf"

echo "ROH analysis completed!"

# Step 5: Summarize total ROH length for pseudogenome
zcat "${savedir}/All_French_ROH_Pseudo_BCFtools" | \
awk -v s="pseudohom" 'BEGIN{sum=0}{if ($2==s){sum+=$6}}END{printf "%s\tROH length: %d\n", s, sum}' > "${savedir}/All_French_Pseudo_ROH_Summary.txt"

zcat "${savedir}/All_Italian_ROH_Pseudo_BCFtools" | \
awk -v s="pseudohom" 'BEGIN{sum=0}{if ($2==s){sum+=$6}}END{printf "%s\tROH length: %d\n", s, sum}' > "${savedir}/All_Italian_Pseudo_ROH_Summary.txt"

echo "ROH summaries saved in:"
echo "${savedir}/All_French_Pseudo_ROH_Summary.txt"
echo "${savedir}/All_Italian_Pseudo_ROH_Summary.txt"
