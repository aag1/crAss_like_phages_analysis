#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)

cd /data/umcg-tifn/crAss_analysis/nt_content

module purge; module load R; module list

Rscript ${dir}/calculate_nt_content.R

chmod 440 CRASS_DB_cl_NT_CONTENT.rds
chmod 440 CRASS_DB_cl_*_SKEW.rds
chmod 440 CRASS_DB_cl_CUMULATIVE_*_SKEW.rds
