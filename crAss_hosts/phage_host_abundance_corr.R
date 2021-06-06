sessionInfo()



### read data
M1 <- read.table(
    '/data/umcg-tifn/crAss_analysis/map_reads/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

M1 <- M1[!grepl('\\|s__', rownames(M1)), ]


M2 <- read.table(
    '/data/umcg-tifn/MetaPhlAn_4cohorts/LLD_LLD2_300OB_IBD_merged_abundance_table.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

colnames(M2) <- sub('_metaphlan$', '', colnames(M2))


t <- read.table(
    '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


z <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/selected_ibd_samples.txt',
    header = FALSE,
    stringsAsFactors = FALSE
)[, 1]



### make calculations
for (k in unique(t$cohort)) {

    # cohort data
    samples <- t$sample_id[ t$cohort == k ]
    if (k == 'IBD') {
        samples <- samples[samples %in% z]
        cat('\n', length(samples), 'IBD samples from subjects w/o stoma were included.\n\n')
    }
    samples <- sapply(samples, function (x) ifelse(grepl('^[0-9]', x), paste0('X', x), x))

    m1 <- M1[, samples]
    m2 <- M2[, samples]


    # retain only taxa found in >10 samples per cohort
    n <- apply(m1, 1, function (v) sum(v > 0))
    m1 <- m1[n > 10, ]

    n <- apply(m2, 1, function (v) sum(v > 0))
    m2 <- m2[n > 10, ]

    N <- nrow(m1) * nrow(m2)


    # correlations
    DF <- data.frame(NULL)

    for (p in rownames(m1)) {

        for (b in rownames(m2)) {

            v1 <- as.numeric(m1[p, ])
            v2 <- as.numeric(m2[b, ])

            OBJ <- cor.test(x = v1, y = v2, method = 'spearman')

            df <- data.frame(
                phage_taxon = p,
                host_taxon = b,
                r = OBJ$estimate,
                p_value = OBJ$p.value,
                FDR = p.adjust(OBJ$p.value, method = 'BH', n = N),
                stringsAsFactors = FALSE
            )
            DF <- rbind(DF, df)

        }

    }

    DF$phage_taxon <- sub('^o__', '', DF$phage_taxon)
    DF$phage_taxon <- gsub('\\|[a-z]__', '|', DF$phage_taxon)

    DF <- DF[order(abs(DF$r), decreasing = TRUE), ]


    # write data
    write.table(
        DF,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE,
        file = paste0(k, '_phage_host_abundance_corr.txt')
    )

}


warnings()
