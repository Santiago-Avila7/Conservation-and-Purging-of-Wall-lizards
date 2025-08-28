#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J BED-polarization
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Polarization/%x_%j.err
#SBATCH -p shared
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT,FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

module load python/3.12.3
module load PDCOLD/23.12
module load pysam/0.22.1-cpeGNU-23.12
 
python3 /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Conservation-and-Purging-of-Wall-lizards/BashScripts/14b_PolarisedSitesBED.py