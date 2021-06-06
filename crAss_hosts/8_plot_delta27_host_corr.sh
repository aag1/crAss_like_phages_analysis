#!/bin/bash
#SBATCH --job-name=job8
#SBATCH --output=job8_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts


module purge; module load R; module list

Rscript ${dir}/plot_delta27_host_corr.R --x_axis_log FALSE

Rscript ${dir}/plot_delta27_host_corr.R --x_axis_log TRUE


chmod 440 delta27_host_corr.pdf
chmod 440 delta27_host_corr_LOG.pdf
