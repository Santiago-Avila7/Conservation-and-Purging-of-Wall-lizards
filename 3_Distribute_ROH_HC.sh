#! /bin/bash -l

# Create a subfolder in the scrip for the temporal jobs
scripdir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/ROH_HC 
jobdir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/ROH_HC 

# Distribute the job per sample and create a job per chromosome.  
while read -r sample; do
    # Create a temporal directory for each sample job
    mkdir -p "${jobdir}/${sample}"
    
    while read -r chr; do
        # Create a copy of the script per chromosome in each sample 
        cp "${scripdir}/3_ROH_HaplotypeCaller.sh" \
           "${jobdir}/${sample}/3_ROH_HaplotypeCaller_${sample}_${chr}.tjob"

        # Replace the sample and chromosome fields in the temporary job
        sed -i "s/sample_field/${sample}/" \
               "${jobdir}/${sample}/3_ROH_HaplotypeCaller_${sample}_${chr}.tjob"

        sed -i "s/chr_field/${chr}/" \
               "${jobdir}/${sample}/3_ROH_HaplotypeCaller_${sample}_${chr}.tjob"

        # Submit the job
        sbatch "${jobdir}/${sample}/3_ROH_HaplotypeCaller_${sample}_${chr}.tjob"
        
    done < Chromosome_names_temp.txt  # Loop through chromosomes
done < Samples_Used_Temp # Loop through samples

# Check jobs:
squeue -u sjaq | nl

# End of script.