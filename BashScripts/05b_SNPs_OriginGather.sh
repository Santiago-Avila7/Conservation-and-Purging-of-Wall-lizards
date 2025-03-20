#!/bin/bash -l

#SBATCH -A naiss2025-22-189
#SBATCH -J MergFilt_SNPs
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/ROH_Genotyping/%x_%j_%a.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/ROH_Genotyping/%x_%j_%a.err
#SBATCH -p shared
#SBATCH -t 18:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Merge the SNPs of all chromosomes and apply filters

#Load modules 
module load bcftools/1.20
module load gatk/4.5.0.0
PICARD_HOME=/pdc/software/eb/software/picard/2.25.5

#Define working directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/ROH_JointGenotyping
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Filtering_steps
filtdir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Filtering_steps
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data
reference=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna

date 

# Merging SNPs of all the chromosomes. 
# French 
java -jar $PICARD_HOME/picard.jar GatherVcfs \
-O ${filtdir}/All_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014743.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014744.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014745.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014746.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014747.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014748.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014749.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014750.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014751.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014752.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014753.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014754.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014755.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014756.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014757.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014758.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014759.1_French_snps_filtered.vcf.gz \
-I ${datadir}/CM014760.1_French_snps_filtered.vcf.gz 

# Index
gatk IndexFeatureFile -I ${filtdir}/All_French_snps_filtered.vcf.gz 
bcftools index ${filtdir}/All_French_snps_filtered.vcf.gz

echo "French Whole-genome  VCF completed" wq

# Italian 
java -jar $PICARD_HOME/picard.jar GatherVcfs \
-O ${filtdir}/All_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014743.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014744.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014745.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014746.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014747.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014748.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014749.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014750.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014751.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014752.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014753.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014754.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014755.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014756.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014757.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014758.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014759.1_Italian_snps_filtered.vcf.gz \
-I ${datadir}/CM014760.1_Italian_snps_filtered.vcf.gz 

# Index
gatk IndexFeatureFile -I ${filtdir}/All_Italian_snps_filtered.vcf.gz 
bcftools index ${filtdir}/All_Italian_snps_filtered.vcf.gz

echo "Italian Whole-genome  VCF completed"

# keep only 'pass' variants and index again (Hard-filtering):
#French 
bcftools view ${filtdir}/All_French_snps_filtered.vcf.gz -f 'PASS,.' -Ob -o ${filtdir}/All_French_snps_filtered_pass.bcf
bcftools index ${filtdir}/All_French_snps_filtered_pass.bcf
#Italian
bcftools view ${filtdir}/All_Italian_snps_filtered.vcf.gz -f 'PASS,.' -Ob -o ${filtdir}/All_Italian_snps_filtered_pass.bcf
bcftools index ${filtdir}/All_Italian_snps_filtered_pass.bcf

echo "Final filters starting"

# Further filtering (Missing genotype data must be specially hard to not bias the Gene diversity analysis)

# For genotype quality (GQ) below 20:
bcftools filter -S . -e 'FMT/GQ<20' ${filtdir}/All_French_snps_filtered_pass.bcf -Ob -o ${filtdir}/All_French_snps_filtered_pass_GQ20.bcf
bcftools index ${filtdir}/All_French_snps_filtered_pass_GQ20.bcf

bcftools filter -S . -e 'FMT/GQ<20' ${filtdir}/All_Italian_snps_filtered_pass.bcf -Ob -o ${filtdir}/All_Italian_snps_filtered_pass_GQ20.bcf
bcftools index ${filtdir}/All_Italian_snps_filtered_pass_GQ20.bcf

echo "Genotype quality applied"

# Remove multialellic sites
bcftools view -M2 ${filtdir}/All_French_snps_filtered_pass_GQ20.bcf -Ob -o ${filtdir}/All_French_snps_filtered_pass_GQ20_biallelic.bcf
bcftools index ${filtdir}/All_French_snps_filtered_pass_GQ20_biallelic.bcf

bcftools view -M2 ${filtdir}/All_Italian_snps_filtered_pass_GQ20.bcf -Ob -o ${filtdir}/All_Italian_snps_filtered_pass_GQ20_biallelic.bcf
bcftools index ${filtdir}/All_Italian_snps_filtered_pass_GQ20_biallelic.bcf

echo "Multiallelic sites removed"

# MAF (5%) - Probably make no sense due to sample numbers
bcftools view -q 0.05:minor ${filtdir}/All_French_snps_filtered_pass_GQ20_biallelic.bcf -Ob -o ${filtdir}/All_French_snps_filtered_pass_GQ20_biallelic_maf.bcf
bcftools index ${filtdir}/All_French_snps_filtered_pass_GQ20_biallelic_maf.bcf

bcftools view -q 0.05:minor ${filtdir}/All_Italian_snps_filtered_pass_GQ20_biallelic.bcf -Ob -o ${filtdir}/All_Italian_snps_filtered_pass_GQ20_biallelic_maf.bcf
bcftools index ${filtdir}/All_Italian_snps_filtered_pass_GQ20_biallelic_maf.bcf

# Delete variants with 70% missing genotype data
bcftools filter -e 'F_MISSING > 0.7' ${filtdir}/All_French_snps_filtered_pass_GQ20_biallelic_maf.bcf -Oz -o ${savedir}/All_French_SNPs.vcf.gz
bcftools stats ${savedir}/All_French_SNPs.vcf.gz > ${savedir}/All_French_SNPs.stats
bcftools --tbi index ${savedir}/All_French_SNPs.vcf.gz

bcftools filter -e 'F_MISSING > 0.7' ${filtdir}/All_Italian_snps_filtered_pass_GQ20_biallelic_maf.bcf -Oz -o ${savedir}/All_Italian_SNPs.vcf.gz
bcftools stats ${savedir}/All_Italian_SNPs.vcf.gz > ${savedir}/All_Italian_SNPs.stats
bcftools --tbi index ${savedir}/All_Italian_SNPs.vcf.gz

echo "Filtering complete" 
date 
