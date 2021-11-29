#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=1gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


module purge; module load R; module list

Rscript plot_abundance_matrices.R

chmod 440 crAss_abundance_matrices.pdf
