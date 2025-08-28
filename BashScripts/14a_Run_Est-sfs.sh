#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Polarization_Pmuralis_Kimmura
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.err
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH --mem=2GB
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

#Environment 
mamba activate est-sfst #Before running 

# Directory 
dir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization"
outdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pmuralis"
mkdir -p ${outdir} 

#Run tests

# Kimmura
est-sfs ${dir}/config-kimura.txt ${dir}/Pmuralis_counts_corr.txt ${dir}/seedfile.txt ${outdir}/PmuralisKimura-Out.txt  ${outdir}/PmuralisKimura-pval.txt