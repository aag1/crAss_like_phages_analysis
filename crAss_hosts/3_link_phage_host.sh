#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts


awk '$13/$14 >= 0.8' blastn_spacersDB1_vs_crAss.out > blastn_spacersDB1_vs_crAss.spacer80match.out
awk '$13/$14 >= 0.8' blastn_spacersDB2_vs_crAss.out > blastn_spacersDB2_vs_crAss.spacer80match.out


module purge; module load R; module list
Rscript ${dir}/link_phage_host.R


chmod 440 blastn_spacersDB*_vs_crAss.spacer80match.out
chmod 440 crAss_host_pairs.txt
