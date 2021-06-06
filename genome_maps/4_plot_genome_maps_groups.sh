#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
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
    --CGC_skewF /data/umcg-tifn/crAss_analysis/nt_content/CRASS_DB_cl_CUMULATIVE_GC_SKEW.rds \
    --AT_skewF /data/umcg-tifn/crAss_analysis/nt_content/CRASS_DB_cl_AT_SKEW.rds \
    --CAT_skewF /data/umcg-tifn/crAss_analysis/nt_content/CRASS_DB_cl_CUMULATIVE_AT_SKEW.rds \
    --groupsF /data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt \
    --outF CRASS_DB_cl_GENOME_MAPS_groups.pdf


chmod 440 CRASS_DB_cl_GENOME_MAPS_groups.pdf
