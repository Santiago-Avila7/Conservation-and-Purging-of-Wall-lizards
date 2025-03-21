#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Ann_SNPs_VCf 
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Purging/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Purging/%x_%j.err
#SBATCH -p shared
#SBATCH -t 6:00:00
#SBATCH --mem=4GB
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Activate the mamaba environment (Be sure that it is active when queueing the job)
mamba activate snpeff
module load bcftools/1.20

# Define directories and genome name
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Origin
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Origin
Genome_Database="PodMur_1.0.99"
Chr_names=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Purging_Chr_names.txt

date 

# Rename chromosome names to match the snpeff database
bcftools annotate --rename-chrs ${Chr_names} ${datadir}/All_French_SNPs.vcf.gz -Oz -o ${savedir}/All_French_Purging_tmp.vcf.gz 
bcftools index --tbi ${savedir}/All_French_Purging_tmp.vcf.gz 

bcftools annotate --rename-chrs ${Chr_names} ${datadir}/All_Italian_SNPs.vcf.gz -Oz -o ${savedir}/All_Italian_Purging_tmp.vcf.gz 
bcftools index --tbi ${savedir}/All_Italian_Purging_tmp.vcf.gz 

# Run SnpEff on the origin subsets
snpEff eff -v ${Genome_Database} ${savedir}/All_French_Purging_tmp.vcf.gz   \
 -stats ${savedir}/All_French_Annotated.html -csvStats ${savedir}/All_French_Annotated.csv > ${savedir}/All_French_annotated.vcf

bgzip ${savedir}/All_French_annotated.vcf
bcftools index --tbi ${savedir}/All_French_Annotated.vcf.gz 

snpEff eff -v ${Genome_Database} ${savedir}/All_Italian_Purging_tmp.vcf.gz   \
 -stats ${savedir}/All_Italian_Annotated.html -csvStats ${savedir}/All_Italian_Annotated.csv > ${savedir}/All_Italian_annotated.vcf

bgzip ${savedir}/All_Italian_annotated.vcf
bcftools index --tbi ${savedir}/All_Italian_Annotated.vcf.gz 

rm ${savedir}/*_Purging_tmp.*

echo "SnpEff annotation complete."

date