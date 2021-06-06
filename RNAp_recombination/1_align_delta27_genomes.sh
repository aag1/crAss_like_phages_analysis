#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/RNAp_recombination/



# -------------------- list all delta27 genomes --------------------
module purge; module load R; module list
Rscript ${dir}/list_all_delta27_genomes.R



# -------------------- get circular delta27 genomes --------------------
base='delta27_genomes'
crass_db_fasta=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB.fasta

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    ${crass_db_fasta} \
    ${base}_all.ids > ${base}_all.fasta

perl /home/umcg-agulyaeva/crAss_analysis/crAss_contigs_detection/identify_circular.pl \
    ${base}_all.fasta \
    ${base}_circular.txt \
    10

cut -d$'\t' -f1 ${base}_circular.txt > ${base}_circular.ids

# Four genomes that are close to completeness and represent vOTUs were included as well.
echo 'cs_ms_40_1' >> ${base}_circular.ids
echo 'NL_crAss001084' >> ${base}_circular.ids
echo 'NL_crAss001078' >> ${base}_circular.ids
echo 'NL_crAss000719' >> ${base}_circular.ids

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    ${crass_db_fasta} \
    ${base}_circular.ids > ${base}_circular.fasta



# -------------------- get last 500 nt of OLNE01000081 --------------------
Rscript ${dir}/get_last_500nt_OLNE01000081.R



# -------------------- BLAST --------------------
module purge; module load BLAST+; module list

makeblastdb \
    -dbtype 'nucl' \
    -in ${base}_circular.fasta \
    -out ${base}_circular

blastn \
    -db ${base}_circular \
    -query OLNE01000081_last_500nt.fasta \
    -out OLNE01000081_last_500nt_VS_circular_delta27_genomes.txt \
    -outfmt 6 \
    -num_threads ${SLURM_CPUS_PER_TASK}



# -------------------- align --------------------
module purge; module load R; module list
Rscript ${dir}/prepare_genomes_for_MSA.R

module purge; module load MAFFT; module list
mafft \
    --thread ${SLURM_CPUS_PER_TASK} \
    ${base}_circular_for_MSA.fasta > ${base}_circular_MSA.fasta



# -------------------- clean up --------------------
rm *ids
rm delta27_genomes_all.fasta
rm delta27_genomes_circular.fasta
rm delta27_genomes_circular.n*
chmod 440 *
