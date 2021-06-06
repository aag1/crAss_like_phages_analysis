#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts



module purge; module load BLAST+; module list
blastn -help > blastn_help.txt

makeblastdb \
    -dbtype 'nucl' \
    -in /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB.fasta \
    -out CRASS_DB

blastn \
    -query /data/umcg-tifn/DATABASES/spacerom_Shmakov_2017/Supplementary_material_spacers.fna \
    -db CRASS_DB \
    -dust no \
    -evalue 1 \
    -task blastn-short \
    -max_target_seqs 1000000 \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen' \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out blastn_spacersDB1_vs_crAss.out



rm CRASS_DB*
chmod 440 blastn_help.txt
chmod 440 blastn_spacersDB1_vs_crAss.out
