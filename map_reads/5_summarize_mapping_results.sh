#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/map_reads


module purge; module load R; module list

for cohort in 'LLD' 'LLD2' '300OB' 'IBD'
do
    Rscript ${dir}/summarize_mapping_results.R \
        --cohort ${cohort} \
        --cl_summary '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl_summary.txt' \
        --read_number '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/'${cohort}'_raw_clean_reads_number.txt' \
        --mapping_dir '/data/umcg-tifn/crAss_analysis/map_reads/map_'${cohort}'_reads'
done


chmod 440 *_crAss_mapping_summary.txt
chmod 440 *_crAss_abundance.txt
chmod 440 *_crAss_depth.rds
