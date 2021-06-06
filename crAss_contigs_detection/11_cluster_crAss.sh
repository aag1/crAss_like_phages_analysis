#!/bin/bash
#SBATCH --job-name=job11
#SBATCH --output=job11_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### output directory
dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/



### cluster
base='CRASS_DB'


module purge; module load MUMmer; module list

nucmer --maxmatch --nooptimize ${base}.fasta ${base}.fasta -t ${SLURM_CPUS_PER_TASK} -p ${base}-nucmer.out

show-coords ${base}-nucmer.out.delta > ${base}-nucmer.out.coords 

export LC_ALL=C
gawk -F ' ' '{if ($12!=$13){print $13" "$0}}' ${base}-nucmer.out.coords | sort > ${base}-nucmer.out.coords.sorted


# script by S. Roux from https://github.com/simroux/ClusterGenomes
module purge; module load Perl; module list

perl /home/umcg-agulyaeva/SOFTWARE/Roux_ClusterGenomes_script/Cluster_genomes_5.1.pl \
    -f ${base}.fasta \
    -i 95 \
    -c 85 \
    -nofna \
    -d 'not_applicable'

rm CRASS_DB-nucmer.out.*
rm CRASS_DB-cover.csv



### table with clustering results
module purge; module load R; module list
Rscript ${dir}/process_cl_results.R



### fasta file with genomes representing clusters
sed '1d' CRASS_DB_cl.txt | cut -d$'\t' -f1 | sort | uniq > ids.txt

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 CRASS_DB.fasta ids.txt > CRASS_DB_cl.fasta

rm ids.txt



### summary table
module purge; module load R; module list

grep '^>' /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/all_genomes.fa | sed 's/^>//' | cut -d' ' -f1 > db673.ids
grep '^>' /data/umcg-tifn/DATABASES/data_Guerin_2018/Data_S1/crass_genomes_249.fasta | sed 's/^>//' | cut -d' ' -f1 > db249.ids
grep '^>' /data/umcg-tifn/DATABASES/crAss_GenBank/crAss-like_NTseq_min10kb_GenBank.fasta | sed 's/^>//' | cut -d' ' -f1 > db146.ids

Rscript ${dir}/build_cl_summary.R

rm *.ids



### permissions
chmod 440 CRASS_DB_95-85.clstr
chmod 440 CRASS_DB_cl*
