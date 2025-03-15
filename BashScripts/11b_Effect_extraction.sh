#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Effect_filter
#SBATCH --output=logs/Purging/%x_%j.out
#SBATCH --error=logs/Purging/%x_%j.err
#SBATCH -p shared
#SBATCH -t 5:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
##SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Use SnpSift to extract informaton about the annotations

conda activate snpeff 
module load bcftools/1.20

# Define directories and genome name
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts

date 

# Filter annotoations based on impact 
SnpSift filter "(ANN[*].IMPACT has 'HIGH')" ${datadir}/All_Origins_annotated.vcf.gz > ${savedir}/All_high.vcf
bgzip ${savedir}/All_high.vcf
bcftools index --tbi ${savedir}/All_high.vcf.gz

SnpSift filter "(ANN[*].IMPACT has 'MODERATE')" ${datadir}/All_Origins_annotated.vcf.gz > ${savedir}/All_moderate.vcf
bgzip ${savedir}/All_moderate.vcf
bcftools index --tbi ${savedir}/All_moderate.vcf.gz

SnpSift filter "(ANN[*].IMPACT has 'lOW')" ${datadir}/All_Origins_annotated.vcf.gz > ${savedir}/All_low.vcf
bgzip ${savedir}/All_low.vcf
bcftools index --tbi ${savedir}/All_low.vcf.gz

SnpSift filter "(ANN[*].IMPACT has 'MODIFIER')" ${datadir}/All_Origins_annotated.vcf.gz > ${savedir}/All_moderate.vcf
bgzip ${savedir}/All_moderate.vcf
bcftools index --tbi ${savedir}/All_moderate.vcf.gz
echo "VCF subsets done"

# Extract Intergenic regions
SnpSift filter "(ANN[0].EFFECT has 'intergenic_region')" ${datadir}/All_Origins_annotated.vcf.gz > ${savedir}/All_intergenic.vcf
bgzip ${savedir}/All_intergenic.vcf
bcftools index --tbi ${savedir}/All_intergenic.vcf.gz
