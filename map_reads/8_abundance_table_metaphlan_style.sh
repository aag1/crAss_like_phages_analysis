#!/bin/bash
#SBATCH --job-name=job8
#SBATCH --output=job8_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/map_reads


module purge; module load R; module list

Rscript ${dir}/abundance_table_metaphlan_style.R


chmod 440 LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt
