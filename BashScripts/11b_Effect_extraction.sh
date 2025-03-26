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
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Origin
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts

date 

# Filter annotoations based on impact for Italian populations
SnpSift filter "(ANN[*].IMPACT has 'HIGH')" ${datadir}/All_Italian_annotated.vcf.gz > ${savedir}/All_Italian_high.vcf
bgzip ${savedir}/All_Italian_high.vcf
bcftools index --tbi ${savedir}/All_Italian_high.vcf.gz

SnpSift filter "(ANN[*].IMPACT has 'MODERATE')" ${datadir}/All_Italian_annotated.vcf.gz > ${savedir}/All_Italian_moderate.vcf
bgzip ${savedir}/All_Italian_moderate.vcf
bcftools index --tbi ${savedir}/All_Italian_moderate.vcf.gz

SnpSift filter "(ANN[*].IMPACT has 'lOW')" ${datadir}/All_Italian_annotated.vcf.gz > ${savedir}/All_Italian_low.vcf
bgzip ${savedir}/All_Italian_low.vcf
bcftools index --tbi ${savedir}/All_Italian_low.vcf.gz

SnpSift filter "(ANN[*].IMPACT has 'MODIFIER')" ${datadir}/All_Italian_annotated.vcf.gz > ${savedir}/All_Italian_moderate.vcf
bgzip ${savedir}/All_Italian_moderate.vcf
bcftools index --tbi ${savedir}/All_moderate.vcf.gz
echo "VCF subsets done"

# Extract Intergenic regions
SnpSift filter "(ANN[0].EFFECT has 'intergenic_region')" ${datadir}/All_Italian_annotated.vcf.gz > ${savedir}/All_Italian_intergenic.vcf
bgzip ${savedir}/All_Italian_intergenic.vcf
bcftools index --tbi ${savedir}/All_Italian_intergenic.vcf.gz
