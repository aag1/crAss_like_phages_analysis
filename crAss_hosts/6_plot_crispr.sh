#!/bin/bash
#SBATCH --job-name=job6
#SBATCH --output=job6_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts


module purge; module load R; module list
Rscript ${dir}/plot_crispr.R


chmod 440 crAss_host_crispr.pdf
