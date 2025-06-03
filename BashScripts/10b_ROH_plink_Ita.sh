#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J ROH_search_Ita
#SBATCH --output=logs/ROHs/%x_%j.out
#SBATCH --error=logs/ROHs/%x_%j.err
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIl
#SBATCH --mail-user=sa5674av-s@student.lu.se

# ROH calling using PLINK - Italian

# Modules needed 
module load bioinfo-tools
module load plink/1.90b4.9
module load bcftools/1.20


# Define input and output directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity/ROHs
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity/ROHs

date 

# Transform to bed - The flag is importance since the chromosome names are not numerical
echo "Changing chromosome names"

bcftools annotate --rename-chrs New_Chr_names.txt ${datadir}/All_Italian_final.vcf -Oz -o ${datadir}/All_Italian_final_tmp.vcf.gz 
bcftools index ${datadir}/All_Italian_final_tmp.vcf.gz

plink --vcf ${datadir}/All_Italian_final_tmp.vcf.gz --make-bed --out ${savedir}/All_Italian_final --allow-extra-chr

echo "Italian bed created"

# Remove temporal files
rm ${datadir}/All_Italian_final_tmp.vcf.gz ${datadir}/All_Italian_final_tmp.vcf.gz.tbi

echo "searching for ROHs"

# Find ROHs 

plink --bfile ${savedir}/All_Italian_final \
      --homozyg \
      --homozyg-window-snp 100 \
      --homozyg-snp 50 \
      --homozyg-kb 500 \
      --homozyg-gap 500 \
      --homozyg-density 60 \
      --homozyg-window-missing 5 \
      --homozyg-window-het 2 \
      --homozyg-window-threshold 0.05\
      --out ${savedir}/All_Italian_ROHs_1

echo "ROH completed"


date 
