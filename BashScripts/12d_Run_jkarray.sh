#!/bin/bash
#SBATCH -A naiss2025-22-189
#SBATCH -J Jk_allele_count9
#SBATCH --array=1-109
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Purging/Jackknife/FrenchvItalian/%x_%j_%a.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/Purging/Jackknife/FrenchvItalian/%x_%j_%a.err
#SBATCH -p shared
#SBATCH -t 3:00:00
#SBATCH --mem=500MB  
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

 
# Load required module
module load python/3.12.3
module load PDC/23.12
module load pysam/0.22.1-cpeGNU-23.12 


## Run the Python script and pass the SLURM array task ID
python3 /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Conservation-and-Purging-of-Wall-lizards/BashScripts/12d_Allele_count_jk.py ${SLURM_ARRAY_TASK_ID}
