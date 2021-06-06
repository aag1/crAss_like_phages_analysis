#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd /data/umcg-tifn/crAss_analysis/build_tree



module purge; module load Anaconda3; module list
source activate IQ-TREE

iqtree \
    -B 1000 \
    -s MSA_crAss_TerL.fasta \
    -T ${SLURM_CPUS_PER_TASK}

conda deactivate



chmod 440 MSA_crAss_TerL.fasta.iqtree
chmod 440 MSA_crAss_TerL.fasta.treefile
chmod 440 MSA_crAss_TerL.fasta.log
