#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Ann_VCf 
#SBATCH --output=logs/Purging/%x_%j.out
#SBATCH --error=logs/Purging/%x_%j.err
#SBATCH -p shared
#SBATCH -t 6:00:00
#SBATCH --mem=4GB
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Annotate SNP dataset using SnpEff

# Activate the mamaba environment (Be sure that it is active when queueing the job)
mamba activate snpeff
module load bcftools/1.20

# Define directories and genome name
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
Genome_Database="PodMur_1.0.99"

date 

# Rename chromosome names to match the snpeff database
bcftools annotate --rename-chrs Purging_Chr_names.txt ${datadir}/All_Origins_SNP_final.vcf -Oz -o ${savedir}/All_Origins_Purging_tmp.vcf.gz 
bcftools index --tbi ${savedir}/All_Origins_Purging_tmp.vcf.gz 

# Run SnpEff on the origin subsets
snpEff eff -v ${Genome_Database} ${savedir}/All_Origins_Purging_tmp.vcf.gz   \
 -stats ${savedir}/All_Origins_Annotated.html -csvStats ${savedir}/All_Origins_Annotated.csv > ${savedir}/All_Origins_annotated.vcf

bgzip ${savedir}/All_Origins_annotated.vcf
bcftools index --tbi ${savedir}/All_Origins_annotated.vcf.gz 

rm ${savedir}/All_Origins_Purging_tmp.*

echo "SnpEff annotation complete."

date