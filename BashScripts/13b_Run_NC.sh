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

module load python/3.12.3
module load PDCOLD/23.12
module load pysam/0.22.1-cpeGNU-23.12

scriptdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts"
inputdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging"
outdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization"

date 

# Verify files exist
required_files=(
    "${scriptdir}/13_Nucliotide_count.py"
    "${inputdir}/All_Origins_annotated.vcf.gz"
    "${scriptdir}/Italian"
    "${scriptdir}/WB_samples"
    "${scriptdir}/WE_samples"
    "${scriptdir}/BU_samples"
)

for file in "${required_files[@]}"; do
    if [[ ! -e "$file" ]]; then
        echo "ERROR: Missing required file $file" >&2
        exit 1
    fi
done

# Run the script
#   - STDOUT → Nucleotide_test.txt
#   - STDERR → SLURM %x_%j.err (no redirection)
python3 "${scriptdir}/13_Nucliotide_count.py" \
    "${inputdir}/All_Origins_annotated.vcf.gz" \
    "${scriptdir}/Italian" \
    "${scriptdir}/WB_samples" \
    "${scriptdir}/WE_samples" \
    "${scriptdir}/BU_samples" \
    > "${outdir}/Nucleotide_test.txt"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Python script failed" >&2
    exit 1
fi

echo "Successfully generated ${outdir}/Nucleotide_test.txt"
date 
