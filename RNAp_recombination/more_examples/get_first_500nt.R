.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()


args <- commandArgs(trailingOnly = TRUE)
genus <- args[1]
genome <- args[2]


L <- read.fasta(
    paste0(genus, '_genomes_all.fasta'),
    seqtype = 'DNA',
    forceDNAtolower = FALSE
)


write.fasta(
    sequences = L[[genome]][ 1:500 ],
    names = paste0(genome, '_first_500nt'),
    nbchar = 80,
    file.out = paste0(genome, '_first_500nt.fasta')
)
