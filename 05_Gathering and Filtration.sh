#!/bin/bash -l

# Merge the SNPs of all chromosomes and apply filters

#Load modules 
module load bcftools/1.20
module load gatk/4.5.0.0
PICARD_HOME=/pdc/software/eb/software/picard/2.25.5

#Define working directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/JointGenotyping

mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs

savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs

#Combination of individual chromosomes VCFs: 

date
echo "Final Gathering starting" 

java -jar $PICARD_HOME/picard.jar GatherVcfs \
-O ${savedir}/All_SNP_filtered.vcf.gz \
-I ${datadir}/CM014743.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014744.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014745.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014746.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014747.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014748.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014749.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014750.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014751.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014752.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014753.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014754.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014755.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014756.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014757.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014758.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014759.1_SNP_filtered.vcf.gz \
-I ${datadir}/CM014760.1_SNP_filtered.vcf.gz 

echo "Whole-genome VCF completed" 

#Index the final VCF 
gatk IndexFeatureFile -I ${savedir}/All_SNP_filtered.vcf.gz
bcftools index ${savedir}/All_SNP_filtered.vcf.gz

# keep only 'pass' SNPs from  hard-filtering and index again:
bcftools view ${savedir}/All_SNP_filtered.vcf.gz -f 'PASS,.' -Ob -o ${savedir}/All_SNP_filtered_pass.bcf
bcftools index ${savedir}/All_SNP_filtered_pass.bcf

echo "Final filters starting"

#further filtering:
#for bi-allelic SNPs:
bcftools view -m2 -M2 -v snps ${savedir}/All_SNP_filtered_pass.bcf -Ob -o ${savedir}/All_SNP_filtered_pass_biallelic.bcf
bcftools index ${savedir}/All_SNP_filtered_pass_biallelic.bcf

#for minor allele frequency:

bcftools view -q 0.03:minor ${savedir}/All_SNP_filtered_pass_biallelic.bcf -Ob -o ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03.bcf
bcftools index ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03.bcf

#for genotype quality (GQ) below 20:
bcftools filter -S . -e 'FMT/GQ<20' ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03.bcf -Ob -o ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03_GQ20.bcf
bcftools index ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03_GQ20.bcf

#delete SNPs within 10 bp of each other:
bcftools +prune -w 10bp -n 1 -N rand ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03_GQ20.bcf -Ob -o ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03_GQ20_thin10.bcf
bcftools index ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03_GQ20_thin10.bcf

#delete SNPs with more than 50% missing genotype data
bcftools filter -e 'F_MISSING > 0.5' ${savedir}/All_SNP_filtered_pass_biallelic_maf0.03_GQ20_thin10.bcf -Ov -o ${savedir}/All_SNP_final.vcf
bcftools stats ${savedir}/All_SNP_final.vcf > ${savedir}/All_SNP_final.stats

echo "Filtering complete" 
date 
