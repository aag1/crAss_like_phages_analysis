.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()



### read data
L <- read.fasta(
    'delta27_genomes_circular.fasta',
    seqtype = 'DNA',
    forceDNAtolower = FALSE
)


tab <- read.table(
    'OLNE01000081_last_500nt_VS_circular_delta27_genomes.txt',
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)
colnames(tab) <- c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')


df <- read.table(
    'delta27_genomes_circular.txt',
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)
endsOverlap <- setNames(df[,2], df[,1])



### process genomes
if (!(all(names(L) %in% c(tab$saccver, 'cs_ms_40_1')))) { stop('A target witout a hit, not cs_ms_40_1!') }

tab <- tab[tab$saccver != 'OLNE01000081', ]

if (nrow(tab) > length(unique(tab$saccver))) { stop('More than one hit to a target!') }


for (i in 1:nrow(tab)) {

    # genome sequence
    ACC <- tab$saccver[i]
    SEQ <- L[[ ACC ]]


    # new start coordinate
    START <- tab$send[i] + 1


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

L[[ 'cs_ms_40_1' ]] <- rev(comp(L[[ 'cs_ms_40_1' ]], forceToLower = FALSE, ambiguous = TRUE))



### write fasta
write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = 'delta27_genomes_circular_for_MSA.fasta'
)
