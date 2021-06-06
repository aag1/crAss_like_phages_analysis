.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(bio3d)
sessionInfo()



### input parameters
option_list = list(
    make_option('--prot'),
    make_option('--tabF'),
	make_option('--aliF'),
    make_option('--seleF')
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
tab <- read.table(
    opt$tabF,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


MSA <- read.fasta(opt$aliF)$ali

if (opt$prot == 'TerL') {
    rownames(MSA) <- sub('_[^_]+_[^_]+_[A-Z]_[0-9]$', '', rownames(MSA))
    rownames(MSA)[ rownames(MSA) == 'crAssphage' ] <- 'NC_024711_crAssphage'
} else {
    rownames(MSA) <- sub('^[^@]+@(.+)_[0-9]+$', '\\1', rownames(MSA))
    rownames(MSA)[ rownames(MSA) == 'alpha_gamma@gene_78_crAssphage' ] <- 'NC_024711_crAssphage'
    rownames(MSA)[ rownames(MSA) == 'beta@gene_16_IAS_virus_KJ003983' ] <- 'IAS_virus_KJ003983'
}


df <- read.table(
    opt$seleF,
    sep = '\t',
    row.names = 1,
    header = FALSE,
    stringsAsFactors = FALSE
)



### aligned proteins of 'old' crAss-like phages
MSA <- MSA[rownames(MSA) %in% tab$cl_repres, ]


non_empty_cols <- which(apply(MSA, 2, function (v) any(v != '-')))
MSA <- MSA[, non_empty_cols]


write.fasta(
    ids = rownames(MSA),
    seqs = MSA,
    file = paste0('old_crAss_', opt$prot, '.fasta')
)



### protein sequenes of 'new' crAss-like phages
library(seqinr)
sessionInfo()


L <- list()


for (i in 1:nrow(tab)) {

    x <- tab$cl_repres[i]
    y <- sub('\\|', '', tab$translation[i])


    if (x %in% rownames(MSA)) { next }
    if (opt$prot == 'portal' & x == 'phage_14_3') { next }


    file1 <- paste0('/data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/', x, '/crAss_key_doms_in_', x, '_', y, '.txt')
    t <- read.table(
        file1,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
    )

    file2 <- paste0('/data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/', x, '/', x, '_', y, '_AA.fasta')
    l <- read.fasta(file2, seqtype = 'DNA', forceDNAtolower = FALSE)


    P <- rownames(t)[ t[, opt$prot] != '-' ]


    if (length(P) == 0) { next }

    if (length(P) == 1) { L[[x]] <- l[[P]] }

    if (length(P) > 1) {

        V <- strsplit(df[x, 1], ',')[[1]]

        S <- c()
        for (z in V) {

            s <- l[[z]]

            if (opt$prot == 'portal' & z == 'ERR844003_ms_1_1') { s <- s[19:length(s)] }
            if (opt$prot == 'portal' & z == 'Woesebacteria_MGFQ01000035_105') { s <- s[21:length(s)] }

            S <- c(S, s)

        }

        L[[x]] <- S

    }

}


write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = paste0('new_crAss_', opt$prot, '.fasta')
)
