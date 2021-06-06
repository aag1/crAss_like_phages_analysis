#!/bin/bash
#SBATCH --job-name=job7
#SBATCH --output=job7_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts


module purge; module load R; module list
Rscript ${dir}/phage_host_abundance_corr.R


chmod 440 *_phage_host_abundance_corr.txt
