.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(IRanges)
sessionInfo()



### read data
t1 <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

t2 <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)



### finalize blastn results
for (db in c('CRASS_DB', 'GenBank')) {

    # read data
    tab <- read.table(
        paste0('OLNE01000081_RNAP_VP02740_orfs_VS_', db, '.txt'),
        sep = '\t',
        header = FALSE,
        stringsAsFactors = FALSE
    )
    colnames(tab) <- c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')


    # add taxonomy
    tab$vOTU <- ''
    tab$genus <- ''
    tab$group <-  ''

    idx <- which(tab$saccver %in% t1$cl_member)

    tab$vOTU[idx] <- sapply(tab$saccver[idx], function (x) t1$cl_repres[t1$cl_member == x])
    tab$genus[idx] <- sapply(tab$vOTU[idx], function (x) t2$genus[t2$cl_repres == x])
    tab$group[idx] <- sapply(tab$vOTU[idx], function (x) t2$group[t2$cl_repres == x])


    # cumulative query coverage by hits to a single target
    tab$qcov_sum <- 0

    for (q in unique(tab$qaccver)) {

        for (s in unique(tab$saccver)) {

            idx <- which((tab$qaccver == q) & (tab$saccver == s))
            if (length(idx) == 0) { next }

            ir <- IRanges::IRanges(start = tab$qstart[idx], end = tab$qend[idx])

            df <- IRanges::reduce(ir)
            df <- as.data.frame(df)
            tab$qcov_sum[idx] <- sum(df$width)

        }

    }


    # write data
    write.table(
        tab,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE,
        file = paste0('OLNE01000081_RNAP_VP02740_orfs_VS_', db, '.final.txt')
    )

}
