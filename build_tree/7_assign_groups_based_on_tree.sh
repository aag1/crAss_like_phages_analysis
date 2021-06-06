#!/bin/bash
#SBATCH --job-name=job7
#SBATCH --output=job7_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/build_tree



module purge; module load R; module list
Rscript ${dir}/assign_groups_based_on_tree.R



chmod 440 MSA_crAss_TerL_MidpointRooted.1.newick
chmod 440 CRASS_DB_cl_summary_TRANSL_GROUPS.1.txt
