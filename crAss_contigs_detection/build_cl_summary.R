sessionInfo()



### read data
dat <- read.table(
        'NL_crAss_contigs_2.txt',
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE)

db673 <- read.table(
        'db673.ids',
        header = FALSE,
        stringsAsFactors = FALSE)[, 1]

db249 <- read.table(
        'db249.ids',
        header = FALSE,
        stringsAsFactors = FALSE)[, 1]

db146 <- read.table(
        'db146.ids',
        header = FALSE,
        stringsAsFactors = FALSE)[, 1]

DF <- read.table(
        'CRASS_DB_cl.txt',
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE)



### build summary table
tab <- data.frame(NULL, stringsAsFactors = FALSE)

for (x in unique(DF$cl_repres)) {

    memb <- DF$cl_member[ DF$cl_repres == x ]

    t <- data.frame(
        cl_repres = x,
        cl_size = length(memb),
        n_LLD   = sum(memb %in% dat$contig_name[dat$cohort_name == 'LLD']),
        n_LLD2  = sum(memb %in% dat$contig_name[dat$cohort_name == 'LLD2']),
        n_300OB = sum(memb %in% dat$contig_name[dat$cohort_name == '300OB']),
        n_IBD   = sum(memb %in% dat$contig_name[dat$cohort_name == 'IBD']),
        n_db673 = sum(memb %in% db673),
        n_db249 = sum(memb %in% db249),
        n_db146 = sum(memb %in% db146),
        stringsAsFactors = FALSE
    )

    tab <- rbind(tab, t)

}

tab <- tab[order(tab$cl_size, decreasing = TRUE), ]



### write summary table
write.table(
    tab,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'CRASS_DB_cl_summary.txt'
)
