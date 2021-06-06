#!/bin/bash
#SBATCH --job-name=job12
#SBATCH --output=job12_%A.out
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### output directory
cd /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/



### build bowtie index
module purge; module load Bowtie2; module list

bowtie2-build \
    CRASS_DB_cl.fasta \
    CRASS_DB_cl \
    --threads ${SLURM_CPUS_PER_TASK}



### generate BED file
module purge; module load EMBOSS; module list

infoseq \
    -auto \
    -nocolumns -noheading \
    -delimiter $'\t0\t' \
    -only -name -length \
    CRASS_DB_cl.fasta > CRASS_DB_cl.bed



### permissions
chmod 440 CRASS_DB_cl*bt2
chmod 440 CRASS_DB_cl.bed
