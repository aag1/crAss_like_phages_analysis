#!/bin/bash


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/RNAp_recombination/more_examples/


module purge; module load R; module list
Rscript ${dir}/recombination_plots.R


chmod 440 RNAp_recombination_alpha_beta_zeta.pdf
