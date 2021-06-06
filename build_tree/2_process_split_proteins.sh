#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/build_tree



### find split proteins
module purge; module load R; module list

Rscript ${dir}/process_split_proteins.R \
    --prot 'TerL' \
    --tabF CRASS_DB_cl_summary_TRANSL.txt \
    --aliF /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/trees/TerL_cMAG.afa

Rscript ${dir}/process_split_proteins.R \
    --prot 'portal' \
    --tabF CRASS_DB_cl_summary_TRANSL.txt \
    --aliF /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/trees/portal.afa



### align split proteins
module purge; module load MAFFT; module list

mafft \
    --add split_crAss_TerL.fasta \
    --thread ${SLURM_CPUS_PER_TASK} \
    repres_crAss_TerL.fasta > MSA_split_crAss_TerL.fasta

mafft \
    --add split_crAss_portal.fasta \
    --thread ${SLURM_CPUS_PER_TASK} \
    repres_crAss_portal.fasta > MSA_split_crAss_portal.fasta



### clean up & permissions
rm repres_crAss_*.fasta
rm split_crAss_*.fasta
chmod 440 split_crAss_*.txt
chmod 440 MSA_split_crAss_*.fasta
