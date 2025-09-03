#!/bin/bash -l

#SBATCH -A naiss2025-22-189
#SBATCH -J PCA_nomiss
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/PopGen/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/PopGen/%x_%j.err
#SBATCH -p shared
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

#PCA using the filtered VCF file. 

#Module 
module load plink/2.00a5.14
module load bcftools/1.20

#Define directories: 
infodir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/PCA
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/PCA

#Change names of chromosomes to valid values: 
bcftools annotate --rename-chrs ${infodir}/New_Chr_names.txt ${datadir}/All_Origins_final.vcf.gz -Oz -o ${savedir}/All_Origins_final_tmp.vcf.gz 
bcftools index --tbi ${savedir}/All_Origins_final_tmp.vcf.gz 

# Convert VCF to PLINK binary format
plink --vcf ${savedir}/All_Origins_final_tmp.vcf.gz  --make-bed --out ${savedir}/All_Origins_final_plink --allow-extra-chr

#Run the PCA
plink -bfile ${savedir}/All_Origins_final_plink --aec --pca --out ${savedir}/All_Origins_PCA_NoMiss

rm ${savedir}/All_Origins_final_tmp.*