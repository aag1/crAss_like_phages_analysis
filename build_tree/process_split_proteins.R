.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(bio3d)
sessionInfo()



### input parameters
option_list = list(
    make_option('--prot'),
    make_option('--tabF'),
	make_option('--aliF')
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



### aligned proteins of representative crAss-like phages
msa <- MSA[c('NC_024711_crAssphage', 'IAS_virus_KJ003983', 'OLOZ01000093', 'OHXK01000005', 'OHFV01000001'), ]


non_empty_cols <- which(apply(msa, 2, function (v) any(v != '-')))
msa <- msa[, non_empty_cols]


write.fasta(
    ids = rownames(msa),
    seqs = msa,
    file = paste0('repres_crAss_', opt$prot, '.fasta')
)



### split proteins of 'new' crAss-like phages
library(seqinr)
sessionInfo()


L <- list()
DF <- data.frame(NULL)


for (i in 1:nrow(tab)) {

    x <- tab$cl_repres[i]
    y <- sub('\\|', '', tab$translation[i])


    if (x %in% rownames(MSA)) { next }


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


    if (length(P) > 1) {

        df <- data.frame(
            V1 = x,
            V2 = paste(P, collapse=','),
            stringsAsFactors = FALSE
        )
        DF <- rbind(DF, df)

        l <- l[P]
        L <- c(L, l)

    }

}


write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    file = paste0('split_crAss_', opt$prot, '.txt')
)


write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = paste0('split_crAss_', opt$prot, '.fasta')
)
