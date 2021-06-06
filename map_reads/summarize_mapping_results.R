.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
sessionInfo()



### input parameters
option_list = list(
    make_option('--cohort'),
	make_option('--cl_summary'),
	make_option('--read_number'),
    make_option('--mapping_dir'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
G <- read.table(
    opt$cl_summary,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE)$cl_repres

N <- read.table(
    opt$read_number,
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)



### summarize data
S <- rownames(N)


DF <- data.frame(NULL)

M <- matrix(
    0,
    nrow = length(G),
    ncol = length(S),
    dimnames = list(G, S)
)

L <- list()


for (x in S) {

    file_cov <- paste0(opt$mapping_dir, '/', x, '/', x, '.coverage.txt')
    file_dep <- paste0(opt$mapping_dir, '/', x, '/', x, '.depth.txt')


    t <- read.table(file_cov, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    idx <- which(t[, 7] > 0.1)


    df <- data.frame(
        sample_id = x,
        num_seq_cov10pct = length(idx),
        pct_reads_mapped = round(sum(t[, 4]) / N[x, 'clean_reads'] * 100, 3),
        pct_reads_mapped_cov10pct = round(sum(t[idx, 4]) / N[x, 'clean_reads'] * 100, 3),
        stringsAsFactors = FALSE
    )
    DF <- rbind(DF, df)


    if (length(idx) == 0) { next }


    M[t[idx, 1], x] <- t[idx, 4] / t[idx, 6] / N[x, 'clean_reads'] * 10^6


    d <- read.table(file_dep, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    for (i in idx) {

        seq_id  <- t[i, 1]
        seq_len <- t[i, 6]

        if (!(seq_id %in% names(L))) { L[[ seq_id ]] <- list() }

        row <- which(d[, 1] == seq_id)
        coo <- d[row, 2]
        cov <- d[row, 3]

        COVERAGE <- rep(0, seq_len)
        COVERAGE[coo] <- cov

        w_center <- seq(from = 501, to = seq_len - 500, by = 200)
        w_from <- w_center - 500
        w_to <- w_center + 500
        w_to[ length(w_to) ] <- seq_len

        L[[ seq_id ]][[ x ]] <- sapply(seq_along(w_center), function (i) {

            mean( COVERAGE[ w_from[i]:w_to[i] ] )

        })

    }

}



### write data
write.table(DF, sep = '\t', quote = FALSE, row.names = FALSE, file = paste0(opt$cohort, '_crAss_mapping_summary.txt'))
write.table(M, sep = '\t', quote = FALSE, file = paste0(opt$cohort, '_crAss_abundance.txt'))
saveRDS(L, file = paste0(opt$cohort, '_crAss_depth.rds'))
