#!/bin/bash



dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/RNAp_recombination/more_examples/



GENUS=('alpha6' 'beta8' 'zeta9')
GENOME=('OLOZ01000098' 'IAS_virus_KJ003983' 'NL_crAss000848')

for i in 0 1 2
do
    genus=${GENUS[$i]}
    genome=${GENOME[$i]}

    # -------------------- list all genus genomes --------------------
    module purge; module load R; module list
    Rscript ${dir}/list_all_genus_genomes.R ${genus}



    # -------------------- get circular genus genomes --------------------
    base1=${genus}'_genomes'
    crass_db_fasta=/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB.fasta

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        ${crass_db_fasta} \
        ${base1}_all.ids > ${base1}_all.fasta

    perl /home/umcg-agulyaeva/crAss_analysis/crAss_contigs_detection/identify_circular.pl \
        ${base1}_all.fasta \
        ${base1}_circular.txt \
        10

    cut -d$'\t' -f1 ${base1}_circular.txt > ${base1}_circular.ids

    # add nearly complete genome
    if [ ${genus} == 'alpha6' ]; then echo 'ERR975093_s_1' >> ${base1}_circular.ids; fi

    /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq \
        -l 80 \
        ${crass_db_fasta} \
        ${base1}_circular.ids > ${base1}_circular.fasta



    # -------------------- get first 500 nt of the selected genome --------------------
    Rscript ${dir}/get_first_500nt.R ${genus} ${genome}



    # -------------------- BLAST --------------------
    base2=${genome}_first_500nt
    module purge; module load BLAST+; module list

    makeblastdb \
        -dbtype 'nucl' \
        -in ${base1}_circular.fasta \
        -out ${base1}_circular

    blastn \
        -db ${base1}_circular \
        -query ${base2}.fasta \
        -out ${base2}_VS_${base1}_circular.txt \
        -outfmt 6 \
        -num_threads 1



    # -------------------- align --------------------
    module purge; module load R; module list
    Rscript ${dir}/prepare_genomes_for_MSA.R ${genus} ${genome}

    # IAS_virus_KJ003983 and NL_crAss001291 vOTUs genomes will be aligned w/o prior reshuffling
    if [ ${genus} == 'beta8' ]; then cp ${base1}_all.fasta ${base1}_circular_for_MSA.fasta; fi

    module purge; module load MAFFT; module list
    mafft \
        --thread 1 \
        ${base1}_circular_for_MSA.fasta > ${base1}_circular_MSA.fasta

done



# -------------------- clean up --------------------
rm *ids
rm *_genomes_all.fasta
rm *_genomes_circular.fasta
rm *_genomes_circular.n*
chmod 440 *
