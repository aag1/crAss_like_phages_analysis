#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/RNAp_recombination/



### get OLNE01000081 orfs RNAP and VP02740
# see /data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/OLNE01000081_TAGq_AA_coordinates.txt

module purge; module load SAMtools; module list

cp -p /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB.fasta .

samtools faidx CRASS_DB.fasta OLNE01000081:50656-66981 > OLNE01000081_RNAP_VP02740_orfs.fasta
sed -i 's/OLNE01000081:50656-66981/OLNE01000081_RNAP/' OLNE01000081_RNAP_VP02740_orfs.fasta

samtools faidx CRASS_DB.fasta OLNE01000081:66984-72821 >> OLNE01000081_RNAP_VP02740_orfs.fasta
sed -i 's/OLNE01000081:66984-72821/OLNE01000081_VP02740/' OLNE01000081_RNAP_VP02740_orfs.fasta



### blastn orfs
module purge; module load BLAST+; module list

makeblastdb \
    -dbtype 'nucl' \
    -in /data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB.fasta \
    -out CRASS_DB

blastn \
    -task 'blastn' \
    -db CRASS_DB \
    -evalue 0.001 \
    -query OLNE01000081_RNAP_VP02740_orfs.fasta \
    -out OLNE01000081_RNAP_VP02740_orfs_VS_CRASS_DB.txt \
    -outfmt 6 \
    -num_threads ${SLURM_CPUS_PER_TASK}

blastn \
    -task 'blastn' \
    -db 'nt' -remote \
    -evalue 0.001 \
    -query OLNE01000081_RNAP_VP02740_orfs.fasta \
    -out OLNE01000081_RNAP_VP02740_orfs_VS_GenBank.txt \
    -outfmt 6



### process blastn results
module purge; module load R; module list

Rscript ${dir}/process_blastn_results.R



### clean up & permissions
rm CRASS_DB.fasta*
rm CRASS_DB.n*
chmod 440 OLNE01000081_RNAP_VP02740_orfs.fasta
chmod 440 OLNE01000081_RNAP_VP02740_orfs_VS_CRASS_DB*
chmod 440 OLNE01000081_RNAP_VP02740_orfs_VS_GenBank*
