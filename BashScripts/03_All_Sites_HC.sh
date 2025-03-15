#!/bin/bash -l
#SBATCH -A naiss2024-22-490
#SBATCH -J ROH_VAR_CALL
#SBATCH --output=%x_%j.out
#SBATCH -p shared
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

#Modules used
module load gatk/4.5.0.0

# Define directories and variables

#Define the sample
sample=sample_field


datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/${sample}/mapping
mkdir -p  /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/${sample}/ROH_HaplotypeCaller
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/${sample}/ROH_HaplotypeCaller

reference=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna

# Define chromosome name (e.g., replace with your desired chromosome name)
chr=chr_field

# Run HaplotypeCaller for the specified chromosome and sample
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R ${reference} \
    -I ${datadir}/${sample}_marked_renamed.bam \
    -L ${chr} \
    -ERC BP_RESOLUTION --output-mode EMIT_ALL_CONFIDENT_SITES \
    -O ${savedir}/${sample}_ROH_${chr}.g.vcf.gz

# Index the output GVCF file
gatk IndexFeatureFile -I ${savedir}/${sample}_ROH_${chr}.g.vcf.gz

date
echo "Variant calling done"