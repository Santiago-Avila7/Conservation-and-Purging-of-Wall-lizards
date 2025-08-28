#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J NucleotideCounts 
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.err
#SBATCH -p shared
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT,FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# It uses a VCF with both ingroup and outgroup samples to get organized nucleotide counts data for est-sfs. 

module load python/3.12.3
module load PDCOLD/23.12
module load pysam/0.22.1-cpeGNU-23.12

scriptdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts"
inputdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pol_VCFs"
outdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization"

date 

# Verify files exist
required_files=(
    "${scriptdir}/13_Nucliotide_count.py"
    "${inputdir}/merged_sites_in_ingroup.vcf.gz"
    "${scriptdir}/Samples_Used"
    "${scriptdir}/Polarization_Out1"
    "${scriptdir}/Polarization_Out2"
    "${scriptdir}/Polarization_Out3"
)

for file in "${required_files[@]}"; do
    if [[ ! -e "$file" ]]; then
        echo "ERROR: Missing required file $file" >&2
        exit 1
    fi
done

# Run the script
python3 "${scriptdir}/13_Nucliotide_count.py" \
    "${inputdir}/merged_sites_in_ingroup.vcf.gz" \
    "${scriptdir}/Samples_Used" \
    "${scriptdir}/Polarization_Out1" \
    "${scriptdir}/Polarization_Out2" \
    "${scriptdir}/Polarization_Out3" \
    > "${outdir}/Pmularis_counts.txt"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Python script failed" >&2
    exit 1
fi

echo "Successfully generated ${outdir}/Pmularis_counts.txt"
date 
