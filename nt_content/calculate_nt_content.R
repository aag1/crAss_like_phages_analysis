.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()



### read data
file <- '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.fasta'
L <- read.fasta(file, forceDNAtolower = FALSE)



### nt content
D <- lapply(names(L), function (s) {

    seq_len <- length(L[[s]])

    w_center <- seq(from = 501, to = seq_len - 500, by = 200)
    w_from <- w_center - 500
    w_to <- w_center + 500
    w_to[ length(w_to) ] <- seq_len
    w_len <- w_to - w_from + 1

    d <- lapply(c('A', 'T', 'G', 'C'), function (b) {

            V <- ifelse(L[[s]] == b, 1, 0)

            sapply(seq_along(w_center), function (i) {

                sum( V[ w_from[i]:w_to[i] ] ) / w_len[i] * 100

            })

    })

    names(d) <- c('A', 'T', 'G', 'C')

    return(d)

})

names(D) <- names(L)

saveRDS(D, file = 'CRASS_DB_cl_NT_CONTENT.rds')



### GC skew
GC_skew <- lapply(D, function (l) {

    ( l[['G']] - l[['C']] ) / ( l[['G']] + l[['C']] )

})

saveRDS(GC_skew, file = 'CRASS_DB_cl_GC_SKEW.rds')



### cumulative GC skew
CGC_skew <- lapply(GC_skew, function (v) {

    idx <- which(!(is.nan(v) | is.infinite(v)))   # for example, FUWD013401386 sequence has long stretch of N letters, resulting in NaN values of GC skew
    w <- rep(NaN, length(v))
    w[idx] <- sapply(1:length(idx), function (i) sum(v[idx][1:i]))
    return(w)

})

saveRDS(CGC_skew, file = 'CRASS_DB_cl_CUMULATIVE_GC_SKEW.rds')



### AT skew
AT_skew <- lapply(D, function (l) {

    ( l[['A']] - l[['T']] ) / ( l[['A']] + l[['T']] )

})

saveRDS(AT_skew, file = 'CRASS_DB_cl_AT_SKEW.rds')



### cumulative AT skew
CAT_skew <- lapply(AT_skew, function (v) {

    idx <- which(!(is.nan(v) | is.infinite(v)))   # for example, FUWD013401386 sequence has long stretch of N letters, resulting in NaN values of AT skew
    w <- rep(NaN, length(v))
    w[idx] <- sapply(1:length(idx), function (i) sum(v[idx][1:i]))
    return(w)

})

saveRDS(CAT_skew, file = 'CRASS_DB_cl_CUMULATIVE_AT_SKEW.rds')
