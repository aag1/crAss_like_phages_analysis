#!/bin/bash
#SBATCH --job-name=job1.1
#SBATCH --output=job1.1_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-1000
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### sample ID
file=/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD_scripts_logs_info/LLD_raw_reads_number_sele.txt
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file} | cut -d$'\t' -f1)
echo '#################### WORKING WITH '${SAMPLE_ID}' SAMPLE ####################'



### output directory
dir1=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection

dir2=/data/umcg-tifn/crAss_analysis/map_reads/map_LLD_reads/${SAMPLE_ID}

if [ -d ${dir2} ]; then rm -rf ${dir2}; fi
mkdir ${dir2}
chmod 750 ${dir2}
cd ${dir2}



### align reads
module purge; module load Bowtie2 SAMtools BEDTools; module list

bowtie2 \
    --very-sensitive \
    -x ${dir1}/CRASS_DB_cl \
    -1 /data/umcg-tifn/LLD_clean_reads/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 /data/umcg-tifn/LLD_clean_reads/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --no-unal --threads ${SLURM_CPUS_PER_TASK} |
samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${SAMPLE_ID}.sorted.bam

samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${SAMPLE_ID}.sorted.bam



### calculate statistics
samtools flagstat ${SAMPLE_ID}.sorted.bam > ${SAMPLE_ID}.flagstat.txt

samtools depth ${SAMPLE_ID}.sorted.bam > ${SAMPLE_ID}.depth.txt

bedtools coverage -a ${dir1}/CRASS_DB_cl.bed -b ${SAMPLE_ID}.sorted.bam > ${SAMPLE_ID}.coverage.txt

chmod 440 *
