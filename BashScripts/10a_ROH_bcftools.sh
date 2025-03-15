#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J ROH_BCFTools
#SBATCH --output=logs/ROHs/%x_%j.out
#SBATCH --error=logs/ROHs/%x_%j.err
#SBATCH -p shared
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIl
#SBATCH --mail-user=sa5674av-s@student.lu.se

module load bcftools/1.20

# Define input and output directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity/ROHs
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity/ROHs

date 

echo "searching for ROHs"

# Find ROHs 

bcftools roh \
 --GTs-only 30 \
 --AF-tag AF \
 -M 1e-8 \
 -Or -o ${savedir}/All_French_ROH_BCFtools \
 ${datadir}/All_French_final.vcf
 
echo "French ROHs found" 

bcftools roh \
 --GTs-only 30 \
 --AF-tag AF \
 -M 1e-8 \
 -Or -o ${savedir}/All_Italian_ROH_BCFtools \
 ${datadir}/All_Italian_final.vcf
 
echo "Italian ROHs found" 

date 
