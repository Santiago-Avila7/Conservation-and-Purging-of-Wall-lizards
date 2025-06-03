#!/bin/bash -l
#SBATCH -A naiss2024-22-490
#SBATCH -J ROH_VAR_CALL
#SBATCH --output=logs/ROH_HC/%x_%j_%a.out
#SBATCH --error=logs/ROH_HC/%x_%j_%a.err
#SBATCH -p shared
#SBATCH -t 8:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=12GB
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se
#SBATCH --array=1-510%40

# Variant calling array (1 job per chromosome for each sample)

# Load GATK module
module load gatk/4.5.0.0

# Define task list file and read the corresponding task
task_file="task_list.txt"
task=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $task_file)

# Extract sample and chromosome from the task
sample=$(echo $task | cut -d ' ' -f 1)
chr=$(echo $task | cut -d ' ' -f 2)

# Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/${sample}/mapping
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/${sample}/ROH_HaplotypeCaller
reference=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna

# Create output directory
mkdir -p $savedir

# Run HaplotypeCaller
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R $reference \
    -I ${datadir}/${sample}_marked_renamed.bam \
    -L $chr \
    -ERC BP_RESOLUTION --output-mode EMIT_ALL_CONFIDENT_SITES \
    -O ${savedir}/${sample}_ROH_${chr}.g.vcf.gz

# Index the output GVCF file
gatk IndexFeatureFile -I ${savedir}/${sample}_ROH_${chr}.g.vcf.gz

date
echo "Variant calling done for sample ${sample}, chromosome ${chr}"
