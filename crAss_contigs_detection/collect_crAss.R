.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(seqinr)
sessionInfo()



### function
build_crAss_table <- function (SAMPLES, FOLDER, DF, COHORT) {

   for (s in SAMPLES) {

        ### table per sample
        tab_file <- paste0(FOLDER, s, '/crAss_in_', s, '.txt')
        if (!file.exists(tab_file)) { next }

        k <- read.table(tab_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)[, 1]

        df <- data.frame(
            cohort_name = COHORT,
            sample_id = s,
            contig_id = k,
            stringsAsFactors = FALSE
        )

        df$length <- as.numeric(lapply(strsplit(df$contig_id, '_'), function (v) v[4]))

        df$ends_overlap <- 0


        ### add ends overlap info
        ends_overlap_file <- paste0(FOLDER, s, '/crAss_in_', s, '_ends_overlap.txt')

        if (file.exists(ends_overlap_file)) {

            q <- read.table(ends_overlap_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)[, 1]

            df$ends_overlap[ df$contig_id %in% q ] <- 1

        }


        ### bind
        DF <- rbind(DF, df)

    }

    return(DF)

}



### build table
f <- c(
    'LLD' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD_scripts_logs_info/LLD_raw_reads_number_sele.txt',
    'LLD2' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/LLD2_scripts_logs_info/LLD2_raw_reads_number.txt',
    '300OB' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/300OB_scripts_logs_info/300OB_raw_reads_number.txt',
    'IBD' = '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/IBD_scripts_logs_info/IBD_raw_reads_number_sele.txt'
)

DF <- data.frame(NULL)
test <- c()

for (i in seq_along(f)) {

    samples <- read.table(f[i], sep = '\t', header = FALSE, stringsAsFactors = FALSE)[, 1]
    test <- c(test, samples)

    cohort <- names(f)[i]

    DF <- build_crAss_table(samples, paste0('crAss_in_', cohort, '/'), DF, cohort)

}


b <- (length(test) == length(unique(test)))
cat('\n\nNo matching sample IDs in different cohorts:', b, '\n\n\n')


DF$contig_name <- sprintf("NL_crAss%06d", 1:nrow(DF))


write.table(
    DF,
    sep = '\t',
    row.names = FALSE,
    quote = FALSE,
    file = 'NL_crAss_contigs.txt'
)



### build fasta
b <- 1

for (s in unique(DF$sample_id)) {

    idx <- which(DF$sample_id == s)
    df <- DF[idx, ]


    seq_file <- paste0('crAss_in_', df$cohort_name[1], '/', s, '/crAss_in_', s, '.fasta')
    L <- read.fasta(seq_file, forceDNAtolower = FALSE)


    names(L) <- sapply(names(L), function (x) { df$contig_name[df$contig_id == x] })


    write.fasta(
        sequences = L,
        names = names(L),
        nbchar = 80,
        file.out = 'NL_crAss_contigs.fasta',
        open = ifelse(b == 1, 'w', 'a')
    )
    b <- 0

}
