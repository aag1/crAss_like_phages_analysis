#!/bin/bash
#SBATCH --job-name=job9
#SBATCH --output=job9_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-520
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### sample ID
file=/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/IBD_scripts_logs_info/IBD_raw_reads_number_sele.txt
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file} | cut -d$'\t' -f1)
echo '#################### WORKING WITH '${SAMPLE_ID}' SAMPLE ####################'



### output directory
dir0=$(pwd)
dir1=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection
dir2=${dir1}/crAss_in_IBD/${SAMPLE_ID}

if [ -d ${dir2} ]; then rm -rf ${dir2}; fi
mkdir ${dir2}
chmod 750 ${dir2}
cd ${dir2}



### get contigs >= 10 kb
module purge; module load BBMap; module list

F1=${SAMPLE_ID}_contigs_min10kb.fasta

reformat.sh \
    in=/data/umcg-tifn/IBD_assembly/${SAMPLE_ID}/${SAMPLE_ID}_contigs.fasta \
    out=${F1} \
    fastawrap=80 \
    minlength=10000

if [ $(stat -c%s ${F1}) -eq 0 ]; then rm ${F1}; fi



######################################## SIX FRAMES ########################################
### translate in 6 frames
module purge; module load EMBOSS; module list

F2a=${SAMPLE_ID}_6frames.fasta

if [ -f ${F1} ]
then
    transeq \
        -frame 6 \
        -clean \
        -sequence ${F1} \
        -outseq ${F2a}
fi



### run hmmer
module purge; module load HMMER/3.3-foss-2018a; module list

F3a=hmmer_crAss_markers_VS_${SAMPLE_ID}_6frames.txt

if [ -f ${F2a} ]
then
    hmmsearch \
        --max \
        -E 0.001 \
        --noali \
        --tblout ${F3a} \
        -o hmmer_out.txt \
        --cpu ${SLURM_CPUS_PER_TASK} \
        ${dir1}/MCP_portal_TerL.hmm \
        ${F2a}

    rm hmmer_out.txt
    
    rm ${F2a}
fi



### parse hmmer output
module purge; module load R; module list

F4a=crAss_in_${SAMPLE_ID}_6frames.txt

if [ -f ${F3a} ] && [ $(grep -c -v '^#' ${F3a}) -gt 0 ]
then
    sed 's/ # / | /g' ${F3a} > ${F3a}_parsable

    Rscript ${dir0}/parse_hmmer_output.R \
        --inF ${F3a}_parsable \
        --outF ${F4a}

    sort -o ${F4a} ${F4a}

    rm ${F3a}_parsable
fi



########################################    ORFs    ########################################
### predict ORFs
module purge; module load prodigal; module list

F2b=${SAMPLE_ID}_orfs.fasta

if [ -f ${F1} ]
then
    prodigal \
        -p meta \
        -i ${F1} \
        -a ${F2b} \
        -o prodigal_out.txt

    rm prodigal_out.txt
fi



### run hmmer
module purge; module load HMMER/3.3-foss-2018a; module list

F3b=hmmer_crAss_markers_VS_${SAMPLE_ID}_orfs.txt

if [ -f ${F2b} ]
then
    hmmsearch \
        --max \
        -E 0.001 \
        --noali \
        --tblout ${F3b} \
        -o hmmer_out.txt \
        --cpu ${SLURM_CPUS_PER_TASK} \
        ${dir1}/MCP_portal_TerL.hmm \
        ${F2b}

    rm hmmer_out.txt
    
    rm ${F2b}
fi



### parse hmmer output
module purge; module load R; module list

F4b=crAss_in_${SAMPLE_ID}_orfs.txt

if [ -f ${F3b} ] && [ $(grep -c -v '^#' ${F3b}) -gt 0 ]
then
    sed 's/ # / | /g' ${F3b} > ${F3b}_parsable

    Rscript ${dir0}/parse_hmmer_output.R \
        --inF ${F3b}_parsable \
        --outF ${F4b}

    sort -o ${F4b} ${F4b}

    rm ${F3b}_parsable
fi



########################################    JOIN    ########################################
### list all identified crAss-like contigs
F5=crAss_in_${SAMPLE_ID}.txt

Qa=0
Qb=0
Q=0

if [ -f ${F4a} ] && [ ! -f ${F4b} ]
then
    cat ${F4a} > ${F5}

    Qa=$(wc -l ${F4a} | cut -d' ' -f1)
    Q=$(wc -l ${F5} | cut -d' ' -f1)
fi

if [ ! -f ${F4a} ] && [ -f ${F4b} ]
then
    cat ${F4b} > ${F5}

    Qb=$(wc -l ${F4b} | cut -d' ' -f1)
    Q=$(wc -l ${F5} | cut -d' ' -f1)
fi

if [ -f ${F4a} ] && [ -f ${F4b} ]
then
    cat ${F4a} > ${F5}
    comm -13 ${F4a} ${F4b} >> ${F5}

    Qa=$(wc -l ${F4a} | cut -d' ' -f1)
    Qb=$(wc -l ${F4b} | cut -d' ' -f1)
    Q=$(wc -l ${F5} | cut -d' ' -f1)
fi

echo -e ${SAMPLE_ID}'\t'${Q}'\t'${Qa}'\t'${Qb} > ${SAMPLE_ID}_crAss_contigs_number.txt



### get crAss-like contigs & identify those with >= 10 nt ends overalap
F6=crAss_in_${SAMPLE_ID}.fasta
F7=crAss_in_${SAMPLE_ID}_ends_overlap.txt

if [ -f ${F5} ]
then
    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq -l 80 ${F1} ${F5} > ${F6}

    perl ${dir0}/identify_circular.pl ${F6} ${F7} 10

    if [ $(stat -c%s ${F7}) -eq 0 ]; then rm ${F7}; fi
fi



### count the number of contigs >= 10 kb
N=0

if [ -f ${F1} ]; then N=$(grep -c '^>' ${F1}); fi

echo -e ${SAMPLE_ID}'\t'${N} > ${SAMPLE_ID}_10kb_contigs_number.txt

rm ${F1}

chmod 440 *
