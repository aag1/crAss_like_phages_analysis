#!/bin/bash
#SBATCH --job-name=job7
#SBATCH --output=job7_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/map_reads


module purge; module load R; module list

Rscript ${dir}/plot_coverage_depth.R \
    --LLD_depth /data/umcg-tifn/crAss_analysis/map_reads/LLD_crAss_depth.rds \
    --LLD2_depth /data/umcg-tifn/crAss_analysis/map_reads/LLD2_crAss_depth.rds \
    --OB_depth /data/umcg-tifn/crAss_analysis/map_reads/300OB_crAss_depth.rds \
    --IBD_depth /data/umcg-tifn/crAss_analysis/map_reads/IBD_crAss_depth.rds \
    --nt_content /data/umcg-tifn/crAss_analysis/nt_content/CRASS_DB_cl_NT_CONTENT.rds \
    --groups_tab /data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt \
    --outF CRASS_DB_cl_DEPTH_groups.pdf


chmod 440 CRASS_DB_cl_DEPTH_groups.pdf
