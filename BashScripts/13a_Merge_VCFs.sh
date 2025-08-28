#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J Merging_VCFs
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.err
#SBATCH -p shared
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT,FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Meging of ingroup and outgroup VCF for polarization. 

# Module
module load bcftools/1.20

# Paths 
INGROUP="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/All_Origins_SNP_final.vcf"
OUTGROUPS_VCF="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Polarization/All_Polarization_final.vcf.gz"
OUTDIR="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pol_VCFs"
REF_FASTA="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna"  

mkdir -p "${OUTDIR}"

# Intermediate files
INGROUP_NOMISS="${OUTDIR}/ingroup.nomiss.vcf.gz"
INGROUP_NORM="${OUTDIR}/ingroup.norm.vcf.gz"        # Normalized ingroup
OUTGROUP_NORM="${OUTDIR}/outgroup.norm.vcf.gz"      # Normalized outgroup
MERGED_ALL="${OUTDIR}/merged_all.vcf.gz"
FINAL="${OUTDIR}/merged_sites_in_ingroup.vcf.gz"

# Step 1: Remove missing genotypes in ingroup
bcftools filter -e 'F_MISSING > 0' "${INGROUP}" -Oz -o "${INGROUP_NOMISS}"
bcftools index -t "${INGROUP_NOMISS}"

# Step 2: Normalize alleles using the reference genome (fixes REF/ALT mismatches)
bcftools norm -f "${REF_FASTA}" -c s "${INGROUP_NOMISS}" -Oz -o "${INGROUP_NORM}"
bcftools norm -f "${REF_FASTA}" -c s "${OUTGROUPS_VCF}" -Oz -o "${OUTGROUP_NORM}"
bcftools index -t "${INGROUP_NORM}"
bcftools index -t "${OUTGROUP_NORM}"

# Step 3: Merge normalized VCFs (preserve all sites)
bcftools merge -m all -Oz -o "${MERGED_ALL}" "${INGROUP_NORM}" "${OUTGROUP_NORM}"
bcftools index -t "${MERGED_ALL}"

# Step 4: Keep only positions from the normalized ingroup
bcftools view -T "${INGROUP_NORM}" "${MERGED_ALL}" -Oz -o "${FINAL}"
bcftools index -t "${FINAL}"

# Step 5: Verify outgroup genotypes aren't missing
bcftools stats "${FINAL}" > "${FINAL%.vcf.gz}.stats"
echo "Check ${FINAL%.vcf.gz}.stats for missingness"

# Cleanup (optional)
# rm -f "${INGROUP_NOMISS}"* "${INGROUP_NORM}"* "${OUTGROUP_NORM}"* "${MERGED_ALL}"*