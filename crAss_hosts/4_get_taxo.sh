#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


dir=$(pwd)
cd /data/umcg-tifn/crAss_analysis/crAss_hosts


arr=( $(sed '1d' crAss_host_pairs.txt | cut -d$'\t' -f2 | sort | uniq) )

module purge; module load Perl/5.26.1-foss-2018a BioPerl/1.6.924-foss-2018a-Perl-5.26.1; module list

perl ${dir}/get_taxo.pl ${arr[@]} > host_taxonomy.txt

chmod 440 host_taxonomy.txt
