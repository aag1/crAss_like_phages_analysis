sessionInfo()



# -------------------- function -----------------
LINK <- function(tab) {

    # new columns
    tab$q_identity <- (tab$nident / tab$qlen) * 100

    tab$s_coo1 <- apply(tab[, c('sstart', 'send')], 1, min)
    tab$s_coo2 <- apply(tab[, c('sstart', 'send')], 1, max)


    # order by bit-score
    tab <- tab[order(tab$bitscore, decreasing = TRUE), ]


    # the same spacer maps to a phage genome more than once?
    t <- tab[, c('qaccver', 'saccver')]
    rownames(t) <- as.character(1:nrow(tab))
    t <- unique(t)
    tab <- tab[as.numeric(rownames(t)), ]


    # multiple spacers of a host map to exactly the same region of a phage genome?
    t <- tab[, c('q_organism', 'saccver', 's_coo1', 's_coo2')]
    rownames(t) <- as.character(1:nrow(tab))
    t <- unique(t)
    tab <- tab[as.numeric(rownames(t)), ]


    # link phages and hosts
    DF <- data.frame(NULL)

    for (phage in unique(tab$saccver)) {

        for (host in unique(tab$q_organism)) {

            idx <- which((tab$saccver == phage) & (tab$q_organism == host))

            idx1 <- idx[ tab$q_identity[idx] >= 95 ]
            idx2 <- idx[ tab$q_identity[idx] >= 80 ]

            SELE <- c()
            if (length(idx1) > 0) { SELE <- idx1 }
            if ((length(idx1) == 0) & (length(idx2) > 1)) { SELE <- idx2 }

            if (length(SELE) > 0) {

                spacer <- paste(tab$qaccver[SELE], collapse = ';')
                df <- data.frame(phage, host, spacer, stringsAsFactors = FALSE)
                DF <- rbind(DF, df)

            }
        }
    }


    # output
    return(DF)

}



# -------------------- Spacer database Shmakov et al. 2017 -----------------
# read data
tab <- read.table(
    'blastn_spacersDB1_vs_crAss.spacer80match.out',
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)

colnames(tab) <- c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'nident', 'qlen')


# process data
tab$q_organism <- unlist(lapply(strsplit(tab$qaccver, '_'), function (v) v[1]))

DF1 <- LINK(tab)



# -------------------- CRISPR-Cas++ spacer database 20210121 -----------------
# read data
t <- read.table(
    'blastn_spacersDB2_vs_crAss.spacer80match.out',
    sep = '\t',
    header = FALSE,
    stringsAsFactors = FALSE
)
colnames(t) <- c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'nident', 'qlen')


d1 <- read.table(
    '/data/umcg-tifn/DATABASES/CRISPRCasdb/20210121_spacer_34_safe_names.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


# process data
d1 <- d1[d1$spacer_id %in% unique(t$qaccver), ]

d2 <- data.frame(NULL)

for (i in 1:nrow(d1)) {

    x <- data.frame(spacer_id = d1$spacer_id[i], q_organism = strsplit(d1$genome_id[i], '\\|')[[1]], stringsAsFactors = FALSE)
    d2 <- rbind(d2, x)

}

tab <- merge(
        x = t,
        y = d2,
        by.x = 'qaccver',
        by.y = 'spacer_id',
        sort = FALSE
)

DF2 <- LINK(tab)



# -------------------- output table -----------------
TAB <- rbind(DF1, DF2)

TAB <- aggregate(TAB$spacer, by = list(TAB$phage, TAB$host), function (v) paste(v, collapse = ';'))
colnames(TAB) <- c('phage', 'host', 'spacer')

write.table(
    TAB,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'crAss_host_pairs.txt'
)
