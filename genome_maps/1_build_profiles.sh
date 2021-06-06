#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


PROT=( $(cut -d$'\t' -f1 selected_domains.txt) )


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/genome_maps


module purge; module load HMMER/3.3-foss-2018a; module list


hmm_file=crAss_key_doms.hmm
if [ -f ${hmm_file} ]; then rm -rf ${hmm_file}; fi


for prot in ${PROT[@]}
do
    arr=( $(awk -v var=${prot} '{if ($2==var) print $1}' /data/umcg-tifn/DATABASES/data_Yutin_2020/profile_list_1.txt) )

    for prof in ${arr[@]}
    do
        name=${prot}'_'${prof}

        if [ ${name} == 'MCP_MCP' ]
        then

            hmmbuild -n ${name} ${name}'.hmm' /data/umcg-tifn/DATABASES/data_Yutin_2017/crassphage_2017/prot_align/MCP.afa

            cat ${name}'.hmm' >> ${hmm_file}

            rm ${name}'.hmm'

        elif [ ${name} == 'BACON_cd14948' ]
        then

            hmmbuild -n ${name} ${name}'.hmm' ${dir}/cd14948.fasta

            cat ${name}'.hmm' >> ${hmm_file}

            rm ${name}'.hmm'

        else

            file1=/data/umcg-tifn/DATABASES/data_Yutin_2020/crassfamily_2020/profiles/${prof}.afac
            file2=${name}'_ids_no_consensus.txt'
            file3=${name}'.fasta'
            file4=${name}'.hmm'

            grep '>' ${file1} | grep -v '_consensus' | sed 's/^>//' > ${file2}
            /home/umcg-agulyaeva/SOFTWARE/seqtk-1.3/seqtk subseq ${file1} ${file2} > ${file3}

            hmmbuild -n ${name} ${file4} ${file3}

            cat ${file4} >> ${hmm_file}

            rm ${file2} ${file3} ${file4}

        fi
    done
done


chmod 440 ${hmm_file}
