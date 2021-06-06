.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()



### read data
tab <- read.table(
    'NL_crAss_contigs_CheckV_contamination.tsv',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

DF <- read.table(
    'NL_crAss_contigs.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$contig_name

L <- read.fasta(
    'NL_crAss_contigs.fasta',
    forceDNAtolower = FALSE
)



### remove contamination recognized by CheckV
tab <- tab[tab$provirus == 'Yes', ]

DF$note <- ''

for (i in 1:nrow(tab)) {

    if (!(tab$contig_id[i] %in% c('NL_crAss000042', 'NL_crAss000541'))) { next }


    if (!(tab$region_types[i] %in% c('viral,host', 'host,viral', 'host,viral,host'))) {

        stop('Unexpected value in \'region_types\' column:', tab$region_types[i], '!')

    }


    S <- tab$contig_id[i]
    reg_type <- strsplit(tab$region_types[i], ',')[[1]]
    reg_coo <- strsplit(tab$region_coords_bp[i], ',')[[1]]


    # remove host regions from the contig sequence
    X <- which(reg_type == 'viral')

    vir_coo <- strsplit(reg_coo[X], '-')[[1]]
    vir_from <- as.numeric(vir_coo[1])
    vir_to <- as.numeric(vir_coo[2])

    vir_seq <- L[[S]][ vir_from:vir_to ]
    L[[S]] <- vir_seq


    # add a note about the removed regions to the table
    Y <- which(reg_type == 'host')

    host_coo <- reg_coo[Y]

    if (length(Y) == 1) {

        DF[S, 'note'] <- paste('Region', host_coo, 'was recognized as contamination by CheckV and discarded.')

    } else {

        DF[S, 'note'] <- paste('Regions', host_coo[1], 'and', host_coo[2], 'were recognized as contamination by CheckV and discarded.')

    }

}



### remove contamination recognized based on read mapping patterns and nucleotide content
S <- 'NL_crAss000622'
vir_from <- 1
vir_to <- 75999
host_coo <- '76000-109148'

vir_seq <- L[[S]][ vir_from:vir_to ]
L[[S]] <- vir_seq

DF[S, 'note'] <- paste('Region', host_coo, 'was recognized as a potential contamination and discarded.')



### write updated table and fasta
write.table(
    DF,
    sep = '\t',
    row.names = FALSE,
    quote = FALSE,
    file = 'NL_crAss_contigs_2.txt'
)

write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = 'NL_crAss_contigs_2.fasta'
)
