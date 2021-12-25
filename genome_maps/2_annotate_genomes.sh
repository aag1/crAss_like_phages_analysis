#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-378
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### genome ID
dir0=$(pwd)
dir1=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection
dir2=/data/umcg-tifn/crAss_analysis/genome_maps

N=$((${SLURM_ARRAY_TASK_ID}+1))
GENOME_ID=$(sed "${N}q;d" ${dir1}/CRASS_DB_cl_summary.txt | cut -d$'\t' -f1)
echo '#################### WORKING WITH '${GENOME_ID}' SAMPLE ####################'



### output directory
dir=${dir2}/ANNOTATIONS/${GENOME_ID}

if [ -d ${dir} ]; then rm -rf ${dir}; fi
mkdir ${dir}
chmod 750 ${dir}
cd ${dir}



### get sequence
echo ${GENOME_ID} > id.txt

/home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
    -l 80 \
    ${dir1}/CRASS_DB_cl.fasta \
    id.txt > ${GENOME_ID}.fasta

rm id.txt



### select prodigal mode
# to use 'single' mode, sequence must be at least 20000 characters
len=$(awk -v id=${GENOME_ID} '{if ($2==id) print $3}' ${dir1}/CRASS_DB_cl.txt)

if [ ${len} -lt 100000 ]
then
    p='meta'
else
    p='single'
fi



### annotate
for m in 'c11' 'TAGq' 'TGAw'
do
    ### predict ORFs
    module purge; module load prodigal; module list

    if [ ${m} == 'c11' ]
    then
        prodigal \
            -p ${p} \
            -i ${GENOME_ID}.fasta \
            -a ${GENOME_ID}_${m}_AA.fasta \
            -o ${GENOME_ID}_${m}_prodigal.out
    fi

    if [ ${m} == 'TAGq' ]
    then
        /home/umcg-agulyaeva/SOFTWARE/prodigal_Yutin_2020/code.prodigal/prodigal \
            -g 11 \
            -TAG Q \
            -p ${p} \
            -i ${GENOME_ID}.fasta \
            -a ${GENOME_ID}_${m}_AA.fasta \
            -o ${GENOME_ID}_${m}_prodigal.out
    fi

    if [ ${m} == 'TGAw' ]
    then
        /home/umcg-agulyaeva/SOFTWARE/prodigal_Yutin_2020/code.prodigal/prodigal \
            -g 11 \
            -TGA W \
            -p ${p} \
            -i ${GENOME_ID}.fasta \
            -a ${GENOME_ID}_${m}_AA.fasta \
            -o ${GENOME_ID}_${m}_prodigal.out
    fi



    ### extract proteome coordinates
    grep '^>' ${GENOME_ID}_${m}_AA.fasta | sed 's/^>//' | sed 's/ \# /\t/g' > ${GENOME_ID}_${m}_AA_coordinates.txt



    ### run hmmer
    module purge; module load HMMER/3.3-foss-2018a; module list

    hmmsearch \
        --max \
        -E 0.001 \
        --cpu ${SLURM_CPUS_PER_TASK} \
        --domtblout ${GENOME_ID}_${m}_hmmer.txt \
        -o ${GENOME_ID}_${m}_hmmer.out \
        ${dir2}/crAss_key_doms.hmm \
        ${GENOME_ID}_${m}_AA.fasta



    ### parse hmmer output
    module purge; module load R; module list

    sed 's/ # / | /g' ${GENOME_ID}_${m}_hmmer.txt > ${GENOME_ID}_${m}_hmmer_parsable.txt

    Rscript ${dir0}/parse_hmmer_output.R \
        --domF ${dir0}/selected_domains.txt \
        --inF ${GENOME_ID}_${m}_hmmer_parsable.txt \
        --outF crAss_key_doms_in_${GENOME_ID}_${m}.txt

    rm ${GENOME_ID}_${m}_hmmer_parsable.txt
done

chmod 440 *
