#! /bin/bash -l
#SBATCH -A naiss2024-22-490
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH -J mapping
#SBATCH -t 15:00:00

date

module load bwa/0.7.17
module load picard/2.25.5
module load samtools/1.20

PICARD_HOME=/pdc/software/eb/software/picard/2.25.5
reference=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna

mkdir fastqcTrim
mkdir mapping

# run fastqc on trimmed files
fastqc -o fastqcTrim/ -f fastq -t 16 trimReads/*

# mapping to the reference genome
bwa mem -M -t 16 ${reference} \
trimReads/Pmu_BNI_P90s_trim_1.fq.gz trimReads/Pmu_BNI_P90s_trim_2.fq.gz > mapping/Pmu_BNI_P90s_pair.sam

bwa mem -M -t 16 ${reference} \
trimReads/Pmu_BNI_P90s_trim_unpair1.fq.gz > mapping/Pmu_BNI_P90s_unpair_1.sam

bwa mem -M -t 16 ${reference} \
trimReads/Pmu_BNI_P90s_trim_unpair2.fq.gz > mapping/Pmu_BNI_P90s_unpair_2.sam

date
echo "bwa done"

# sam to bam
java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups -INPUT mapping/Pmu_BNI_P90s_pair.sam -OUTPUT mapping/Pmu_BNI_P90s_pair.bam -SORT_ORDER coordinate \
-RGLB 1708 -RGPL ILLUMINA -RGPU 386 -RGSM Pmu_BNI_P90s

java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups -INPUT mapping/Pmu_BNI_P90s_unpair_1.sam -OUTPUT mapping/Pmu_BNI_P90s_unpair_1.bam -SORT_ORDER coordinate \
-RGLB 1708 -RGPL ILLUMINA -RGPU 386 -RGSM Pmu_BNI_P90s

java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups -INPUT mapping/Pmu_BNI_P90s_unpair_2.sam -OUTPUT mapping/Pmu_BNI_P90s_unpair_2.bam -SORT_ORDER coordinate \
-RGLB 1708 -RGPL ILLUMINA -RGPU 386 -RGSM Pmu_BNI_P90s

# merge bam
java -jar $PICARD_HOME/picard.jar MergeSamFiles -INPUT mapping/Pmu_BNI_P90s_pair.bam -INPUT mapping/Pmu_BNI_P90s_unpair_1.bam -INPUT mapping/Pmu_BNI_P90s_unpair_2.bam \
-OUTPUT mapping/Pmu_BNI_P90s.bam -SORT_ORDER coordinate -ASSUME_SORTED true

# index bam
java -jar $PICARD_HOME/picard.jar BuildBamIndex -INPUT mapping/Pmu_BNI_P90s.bam

# remove intermediate files
rm mapping/Pmu_BNI_P90s_*

#check the mapping stats:
samtools flagstat mapping/Pmu_BNI_P90s.bam > mapping/Pmu_BNI_P90s_MappingStats

date
echo "intermediate done"

#mark duplicates
java -jar $PICARD_HOME/picard.jar MarkDuplicates -INPUT mapping/Pmu_BNI_P90s.bam -OUTPUT mapping/Pmu_BNI_P90s_marked.bam -METRICS_FILE mapping/Pmu_BNI_P90s_marked.metrics

#Index .bam file again
java -jar $PICARD_HOME/picard.jar BuildBamIndex -INPUT mapping/Pmu_BNI_P90s_marked.bam

# remove intermediate files
rm mapping/Pmu_BNI_P90s.ba*

#add meaningful IDs
java -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups -INPUT mapping/Pmu_BNI_P90s_marked.bam -OUTPUT mapping/Pmu_BNI_P90s_marked_renamed.bam -RGID Pmu_BNI_P90s -RGLB Pmu_BNI_P90s-lib \
-RGPL ILLUMINA -RGPU Pmu_BNI_P90s-01 -RGSM Pmu_BNI_P90s

#index again
java -jar $PICARD_HOME/picard.jar BuildBamIndex -INPUT mapping/Pmu_BNI_P90s_marked_renamed.bam

# check the coverage stats
samtools coverage mappingPmu_BNI_P90s_marked_renamed.bam  grep -v '#'  awk '{sum+=$7} END { print ${ind}, AverageCoverage = ,sumNR}'  mappingPmu_BNI_P90s_CoverageStats

# remove intermediate files
rm mapping/Pmu_BNI_P90s_marked.ba*

date
echo "all done and fine"
