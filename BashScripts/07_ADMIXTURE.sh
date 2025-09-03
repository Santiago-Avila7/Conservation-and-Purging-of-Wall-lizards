#!/bin/bash -l

#SBATCH -A naiss2025-22-189
#SBATCH -J ADMIXTURE 
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/PopGen/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/PopGen/%x_%j.err
#SBATCH -p shared
# Mem can be 32GB for K1 to K4
#SBATCH --mem=64GB
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=sa5674av-s@student.lu.se

# Evaluate population structure for K 1 to 10 

# Modules
module load  bioinfo-tools
module load ADMIXTURE/1.3.0

#Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Admixture 
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Admixture 

#Run ADMIXTURE for 1 to 10 ks
for K in {1..10}; do
    echo "Running ADMIXTURE for K=${K}"
    admixture --cv -j8 ${datadir}/All_Origins_final_plink.bed $K | tee ${savedir}/Admixture_K${K}.out
done

# Gather CV errors into a summary file
grep -h "CV error" ${savedir}/Admixture_K*.out > ${savedir}/CV_errors_summary.txt
echo "Cross-validation completed. Results saved in ${savedir}/CV_errors_summary.txt"

#End 
