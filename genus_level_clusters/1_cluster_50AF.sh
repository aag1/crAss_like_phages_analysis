#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### output directory
cd /data/umcg-tifn/crAss_analysis/genus_level_clusters/



### cluster
cp /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.fasta .

base='CRASS_DB_cl'


module purge; module load MUMmer; module list

nucmer --maxmatch --nooptimize ${base}.fasta ${base}.fasta -t ${SLURM_CPUS_PER_TASK} -p ${base}-nucmer.out

show-coords ${base}-nucmer.out.delta > ${base}-nucmer.out.coords 

export LC_ALL=C
gawk -F ' ' '{if ($12!=$13){print $13" "$0}}' ${base}-nucmer.out.coords | sort > ${base}-nucmer.out.coords.sorted


# script by S. Roux from https://github.com/simroux/ClusterGenomes
module purge; module load Perl; module list

perl /home/umcg-agulyaeva/SOFTWARE/Roux_ClusterGenomes_script/Cluster_genomes_5.1.pl \
    -f ${base}.fasta \
    -i 0 \
    -c 50 \
    -nofna \
    -d 'not_applicable'



### clean up & permissions
rm ${base}.fasta
rm ${base}-nucmer.out.*
rm ${base}-cover.csv
chmod 440 ${base}_0-50.clstr
