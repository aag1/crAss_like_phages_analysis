#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/RNAp_recombination/



### Annotate OLNE01000081 proteome
module purge; module load HMMER/3.3-foss-2018a; module list


DB='crassfamily_profiles.hmm'

for file1 in /data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/profiles/*afac
do
    base=$(basename ${file1} .afac)
    file2=${base}'_ids_no_consensus.txt'
    file3=${base}'.fasta'
    file4=${base}'.hmm'

    grep '>' ${file1} | grep -v '_consensus' | sed 's/^>//' > ${file2}
    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq ${file1} ${file2} > ${file3}

    hmmbuild -n ${base} ${file4} ${file3}

    cat ${file4} >> ${DB}

    rm ${file2} ${file3} ${file4}
done


base_out='OLNE01000081_VS_crassfamily_profiles'

hmmsearch \
    --max \
    -E 0.001 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    --domtblout ${base_out}.txt \
    -o ${base_out}.out \
    ${DB} \
    /data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/OLNE01000081_TAGq_AA.fasta

head -n3 ${base_out}.txt > ${base_out}.SELE.txt
grep -E 'OLNE01000081_27|OLNE01000081_28|OLNE01000081_40' ${base_out}.txt >> ${base_out}.SELE.txt



### Find IAS homologs of OLNE01000081 proteins
module purge; module load BLAST+; module list


makeblastdb \
    -dbtype 'prot' \
    -in /data/umcg-tifn/DATABASES/data_Yutin_2017/crassphage_2017/IAS_KJ003983.fa \
    -out IAS_proteome


blastp \
    -db IAS_proteome \
    -query /data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/OLNE01000081_TAGq_AA.fasta \
    -out OLNE01000081_VS_IAS_proteome.txt \
    -outfmt 6 \
    -num_threads ${SLURM_CPUS_PER_TASK}

grep -E 'OLNE01000081_27|OLNE01000081_28|OLNE01000081_40' OLNE01000081_VS_IAS_proteome.txt > OLNE01000081_VS_IAS_proteome.SELE.txt



### Plot
module purge; module load R; module list
Rscript ${dir}/recombination_plot.R



### Clean up & permissions
rm crassfamily_profiles.hmm
rm OLNE01000081_TAGq_VS_crassfamily_profiles.out
rm IAS_proteome.p*
chmod 440 OLNE01000081_VS_crassfamily_profiles*
chmod 440 OLNE01000081_VS_IAS_proteome*
chmod 440 delta27_vOTUs_palette.rds
chmod 440 RNAp_recombination_plot.pdf
