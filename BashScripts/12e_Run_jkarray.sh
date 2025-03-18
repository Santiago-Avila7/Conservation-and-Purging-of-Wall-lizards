#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Jk_allele_count
#SBATCH --array=1-109
#SBATCH --output=logs/Purging/Jackknife/%x_%j_%a.out
#SBATCH --error=logs/Purging/Jackknife/%x_%j_%a.err
#SBATCH -p shared
#SBATCH -t 5:00:00
#SBATCH --mem=1GB  
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Load required module
module load python/3.12.3

## Run the Python script and pass the SLURM array task ID
python3 /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/12e_New_jkarray.py ${SLURM_ARRAY_TASK_ID}