#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J VCF-Filtering
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.err
#SBATCH -p shared
#SBATCH -t 2:00:00  
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT,FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Produced a filtered VCF with only ancestral REF alleles for annotation and downstream analyses 

# Load modules
module load bcftools/1.20

# Set variables
INPUT_VCF="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/All_Origins_SNP_final.vcf"
BED_FILE="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pmuralis/high_confidence_ref_ancestral.bed"
OUTPUT_VCF="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pmuralis/All_Origins_ann_pol.vcf.gz"

echo "Starting VCF filtering"

# Filter VCF using bcftools
bcftools view -T ${BED_FILE} ${INPUT_VCF} -Oz -o ${OUTPUT_VCF}

# Index the filtered VCF
bcftools index --tbi ${OUTPUT_VCF}

echo "VCF filtering complete!"
echo "Filtered VCF created: ${OUTPUT_VCF}"
