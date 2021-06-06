#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### output directory
dir0=$(pwd)
dir1=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection
dir2=${dir1}/crAss_in_db146

if [ -d ${dir2} ]; then rm -rf ${dir2}; fi
mkdir ${dir2}
chmod 750 ${dir2}
cd ${dir2}



### genomes
ntF=/data/umcg-tifn/DATABASES/crAss_GenBank/crAss-like_NTseq_min10kb_GenBank.fasta



######################################## SIX FRAMES ########################################
### translate in 6 frames
module purge; module load EMBOSS; module list

aaF=db146_6frames.fasta

transeq \
    -frame 6 \
    -clean \
    -sequence ${ntF} \
    -outseq ${aaF}



### run hmmer
module purge; module load HMMER/3.3-foss-2018a; module list

tabF=hmmer_crAss_markers_VS_db146_6frames.txt

hmmsearch \
    --max \
    -E 0.001 \
    --noali \
    --tblout ${tabF} \
    -o hmmer_out.txt \
    --cpu ${SLURM_CPUS_PER_TASK} \
    ${dir1}/MCP_portal_TerL.hmm \
    ${aaF}

rm hmmer_out.txt

rm ${aaF}



### parse hmmer output
module purge; module load R; module list

tmpF=hmmer_crAss_markers_VS_db146_6frames_PARSABLE.txt
sed 's/ # / | /g' ${tabF} > ${tmpF}

outF_6fr=crAss_in_db146_6frames.txt

Rscript ${dir0}/parse_hmmer_output.R \
    --inF ${tmpF} \
    --outF ${outF_6fr}

sort -o ${outF_6fr} ${outF_6fr}

rm ${tmpF}



########################################    ORFs    ########################################
### predict ORFs
module purge; module load prodigal; module list

aaF=db146_orfs.fasta

prodigal \
    -p meta \
    -i ${ntF} \
    -a ${aaF} \
    -o prodigal_out.txt

rm prodigal_out.txt



### run hmmer
module purge; module load HMMER; module list

tabF=hmmer_crAss_markers_VS_db146_orfs.txt

hmmsearch \
    --max \
    -E 0.001 \
    --noali \
    --tblout ${tabF} \
    -o hmmer_out.txt \
    --cpu ${SLURM_CPUS_PER_TASK} \
    ${dir1}/MCP_portal_TerL.hmm \
    ${aaF}

rm hmmer_out.txt

rm ${aaF}



### parse hmmer output
module purge; module load R; module list

tmpF=hmmer_crAss_markers_VS_db146_orfs_PARSABLE.txt
sed 's/ # / | /g' ${tabF} > ${tmpF}

outF_orfs=crAss_in_db146_orfs.txt

Rscript ${dir0}/parse_hmmer_output.R \
    --inF ${tmpF} \
    --outF ${outF_orfs}

sort -o ${outF_orfs} ${outF_orfs}

rm ${tmpF}



########################################    JOIN    ########################################
cat ${outF_6fr} > crAss_in_db146.txt
comm -13 ${outF_6fr} ${outF_orfs} >> crAss_in_db146.txt



### database entries that were NOT detected
sort -o crAss_in_db146.txt crAss_in_db146.txt

grep '>' ${ntF} | sed 's/^>//' | cut -d' ' -f1 | sort > db146_all_ids.txt

comm -23 db146_all_ids.txt crAss_in_db146.txt > db146_undetected_ids.txt

rm db146_all_ids.txt



### identify entries with ends overlap
perl ${dir0}/identify_circular.pl ${ntF} db146_ends_overlap.txt 10



### permissions
chmod 440 *
