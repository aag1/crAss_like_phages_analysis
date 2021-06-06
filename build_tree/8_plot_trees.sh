#!/bin/bash
#SBATCH --job-name=job8
#SBATCH --output=job8_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/build_tree



module purge; module load R; module list
Rscript ${dir}/plot_trees.R



chmod 440 MSA_crAss_TerL_MidpointRooted.2.newick
chmod 440 CRASS_DB_cl_summary_TRANSL_GROUPS.2.txt
chmod 440 CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt
chmod 440 crAss_tree_TerL.pdf
chmod 440 MSA_crAss_portal_MidpointRooted.newick
chmod 440 crAss_tree_portal.pdf
