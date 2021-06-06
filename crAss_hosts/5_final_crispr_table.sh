#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts


module purge; module load R; module list
Rscript ${dir}/final_crispr_table.R


chmod 440 crAss_host_pairs.final.txt
