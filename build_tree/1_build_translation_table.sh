#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/build_tree


module purge; module load R; module list
Rscript ${dir}/build_translation_table.R


chmod 440 CRASS_DB_cl_summary_TRANSL.txt
