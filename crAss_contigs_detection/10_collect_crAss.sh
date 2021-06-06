#!/bin/bash
#SBATCH --job-name=job10
#SBATCH --output=job10_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### output directory
dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/



### collect info
echo -e 'sample_id\tn_contigs' > 10kb_contigs_number.txt
cat /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/crAss_in_*/*/*_10kb_contigs_number.txt >> 10kb_contigs_number.txt

echo -e 'sample_id\tn_crAss\tn_crAss_6frames\tn_crAss_orfs' > crAss_contigs_number.txt
cat /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/crAss_in_*/*/*_crAss_contigs_number.txt >> crAss_contigs_number.txt



### build table and fasta of crAss-like contigs
module purge; module load R; module list
Rscript ${dir}/collect_crAss.R



### remove potential host fragments from crAss-like contigs
module purge; module load Anaconda3; module list

source activate CheckV
conda list

checkv contamination \
    NL_crAss_contigs.fasta \
    NL_crAss_contigs_CheckV \
    -d /data/umcg-tifn/DATABASES/checkv-db-v0.6 \
    -t ${SLURM_CPUS_PER_TASK}

conda deactivate


mv NL_crAss_contigs_CheckV/contamination.tsv NL_crAss_contigs_CheckV_contamination.tsv
rm -rf NL_crAss_contigs_CheckV


module purge; module load R; module list
Rscript ${dir}/remove_contamination.R



### combine databases
cat NL_crAss_contigs_2.fasta > CRASS_DB.fasta

db673=/data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/all_genomes.fa
db249=/data/umcg-tifn/DATABASES/data_Guerin_2018/Data_S1/crass_genomes_249.fasta
db146=/data/umcg-tifn/DATABASES/crAss_GenBank/crAss-like_NTseq_min10kb_GenBank.fasta

for db1 in ${db673} ${db249} ${db146}
do
    grep '^>' ${db1} | sed 's/^>//' | cut -d' ' -f1 | sort > db1.ids

    grep '^>' CRASS_DB.fasta | sed 's/^>//' | cut -d' ' -f1 | sort > db2.ids

    comm -23 db1.ids db2.ids > db1_unique.ids

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${db1} db1_unique.ids >> CRASS_DB.fasta

    rm *.ids
done



### permissions
chmod 440 NL_crAss_contigs*
chmod 440 10kb_contigs_number.txt
chmod 440 crAss_contigs_number.txt
chmod 440 CRASS_DB.fasta
