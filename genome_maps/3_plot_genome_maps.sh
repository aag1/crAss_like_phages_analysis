#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/genome_maps


module purge; module load R; module list

Rscript ${dir}/plot_genome_maps.R \
    --domainsF ${dir}/selected_domains.txt \
    --genomesF /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt \
    --GC_skewF /data/umcg-tifn/crAss_analysis/nt_content/CRASS_DB_cl_GC_SKEW.rds \
    --cGC_skewF /data/umcg-tifn/crAss_analysis/nt_content/CRASS_DB_cl_CUMULATIVE_GC_SKEW.rds \
    --outF CRASS_DB_cl_GENOME_MAPS.pdf


chmod 440 CRASS_DB_cl_GENOME_MAPS.pdf
