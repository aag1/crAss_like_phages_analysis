#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/build_tree



module purge; module load R; module list
Rscript ${dir}/prepare_to_align.R \
    --prot 'TerL' \
    --tabF CRASS_DB_cl_summary_TRANSL.txt \
    --aliF /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/trees/TerL_cMAG.afa \
    --seleF /home/umcg-agulyaeva/crAss_analysis/build_tree/SELE_split_crAss_TerL.txt



module purge; module load MAFFT; module list
mafft \
    --add new_crAss_TerL.fasta \
    --thread ${SLURM_CPUS_PER_TASK} \
    old_crAss_TerL.fasta > MSA_crAss_TerL.fasta



chmod 440 old_crAss_TerL.fasta
chmod 440 new_crAss_TerL.fasta
chmod 440 MSA_crAss_TerL.fasta
