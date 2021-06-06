.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()



args <- commandArgs(trailingOnly = TRUE)
genus <- args[1]
genome <- args[2]



### read data
L <- read.fasta(
    paste0(genus, '_genomes_circular.fasta'),
    seqtype = 'DNA',
    forceDNAtolower = FALSE
)


tab <- read.table(
    paste0(genome, '_first_500nt_VS_', genus, '_genomes_circular.txt'),
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)
colnames(tab) <- c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')


df <- read.table(
    paste0(genus, '_genomes_circular.txt'),
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)
endsOverlap <- setNames(df[,2], df[,1])



### process genomes
tab <- tab[tab$saccver != genome, ]

if (genus == 'alpha6') { tab <- tab[-c(7, 10), ] }

if (nrow(tab) > length(unique(tab$saccver))) { stop('More than one hit to a target!') }


for (i in 1:nrow(tab)) {

    # genome sequence
    ACC <- tab$saccver[i]
    SEQ <- L[[ ACC ]]


    # new start coordinate
    START <- tab$sstart[i] + 1


    # reverse complement sequence
    if (tab$sstart[i] > tab$send[i]) {

        SEQ <- rev(comp(SEQ, forceToLower = FALSE, ambiguous = TRUE))

        START <- length(SEQ) - START + 1
    }


    # remove 3'-end nucleotides matching the 5'-end
    if (ACC %in% names(endsOverlap)) {

        SEQ <- SEQ[ 1:(length(SEQ) - endsOverlap[ACC]) ]

        if (length(SEQ) < START) { stop('After removing the 3\'-end nucleotides matching the 5\'-end, the new start coordinate is gone!') }
    }


    # shift the start coordinate
    SEQ <- SEQ[ c(START:length(SEQ), 1:(START-1)) ]


    # update genome
    L[[ ACC ]] <- SEQ

}



### write fasta
write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = paste0(genus, '_genomes_circular_for_MSA.fasta')
)
