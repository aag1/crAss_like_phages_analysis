#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### output directory
dir0=$(pwd)
dir1=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection
dir2=${dir1}/crAss_in_viral_refseq

if [ -d ${dir2} ]; then rm -rf ${dir2}; fi
mkdir ${dir2}
chmod 750 ${dir2}
cd ${dir2}



### genomes
ntF=/data/umcg-tifn/DATABASES/viral_refseq_202/viral_refseq_202_genomic.fna



######################################## SIX FRAMES ########################################
### translate in 6 frames
module purge; module load EMBOSS; module list

aaF=viral_refseq_202_genomic_6frames.fasta

transeq \
    -frame 6 \
    -clean \
    -sequence ${ntF} \
    -outseq ${aaF}



### run hmmer
module purge; module load HMMER/3.3-foss-2018a; module list

tabF=hmmer_crAss_markers_VS_viral_refseq_6frames.txt

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

tmpF=hmmer_crAss_markers_VS_viral_refseq_6frames_PARSABLE.txt
sed 's/ # / | /g' ${tabF} > ${tmpF}

outF_6fr=crAss_in_viral_refseq_6frames.txt

Rscript ${dir0}/parse_hmmer_output.R \
    --inF ${tmpF} \
    --outF ${outF_6fr}

sort -o ${outF_6fr} ${outF_6fr}

rm ${tmpF}



########################################    ORFs    ########################################
### predict ORFs
module purge; module load prodigal; module list

aaF=viral_refseq_202_genomic_orfs.fasta

prodigal \
    -p meta \
    -i ${ntF} \
    -a ${aaF} \
    -o prodigal_out.txt

rm prodigal_out.txt



### run hmmer
module purge; module load HMMER; module list

tabF=hmmer_crAss_markers_VS_viral_refseq_orfs.txt

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

tmpF=hmmer_crAss_markers_VS_viral_refseq_orfs_PARSABLE.txt
sed 's/ # / | /g' ${tabF} > ${tmpF}

outF_orfs=crAss_in_viral_refseq_orfs.txt

Rscript ${dir0}/parse_hmmer_output.R \
    --inF ${tmpF} \
    --outF ${outF_orfs}

sort -o ${outF_orfs} ${outF_orfs}

rm ${tmpF}



########################################    JOIN    ########################################
cat ${outF_6fr} > crAss_in_viral_refseq.txt
comm -13 ${outF_6fr} ${outF_orfs} >> crAss_in_viral_refseq.txt



### permissions
chmod 440 *
