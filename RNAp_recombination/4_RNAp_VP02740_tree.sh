#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/RNAp_recombination/



### get RNAP and VP02740 MSA fragment
module purge; module load R; module list
Rscript ${dir}/get_MSA_fragment.R



### build tree
module purge; module load Anaconda3; module list
source activate IQ-TREE

iqtree \
    -B 1000 \
    -s delta27_RNAP_VP02740_MSA.fasta \
    -T ${SLURM_CPUS_PER_TASK}

conda deactivate



### plot tree
module purge; module load R; module list
Rscript ${dir}/plot_tree.R



rm *.uniqueseq.phy *.model.gz *.mldist *.bionj *splits.nex *.contree *.ckp.gz
chmod 440 delta27_RNAP_VP02740_MSA.fasta*
chmod 440 delta27_RNAP_VP02740_tree.pdf
