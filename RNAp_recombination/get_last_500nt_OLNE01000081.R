.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()


L <- read.fasta(
    'delta27_genomes_all.fasta',
    seqtype = 'DNA',
    forceDNAtolower = FALSE
)

len <- length(L[['OLNE01000081']])


write.fasta(
    sequences = L[['OLNE01000081']][ (len-499):len ],
    names = 'OLNE01000081_last_500nt',
    nbchar = 80,
    file.out = 'OLNE01000081_last_500nt.fasta'
)
