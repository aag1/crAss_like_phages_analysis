sessionInfo()



### read data
h <- read.table(
    'crAss_host_pairs.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


t <- read.table(
    'host_taxonomy.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


k <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt',
    sep = '\t',
    header = TRUE,
    row.names = 2,
    stringsAsFactors = FALSE
)


g <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)



### merge data
df <- merge(
        x = h[, c('phage', 'host')],
        y = t[, c('genome_id', 'taxonomy')],
        by.x = 'host',
        by.y = 'genome_id',
        all.x = TRUE,
        all.y = TRUE,
        sort = FALSE
)
colnames(df) <- c('host_id', 'phage_id', 'host_taxonomy')


df$phage_vOTU <- k[df$phage_id, 'cl_repres']


df$phage_genus <- g[df$phage_vOTU, 'genus']
df$phage_group <- g[df$phage_vOTU, 'group']


df <- df[, c('phage_id', 'phage_vOTU', 'phage_genus', 'phage_group', 'host_id', 'host_taxonomy')]


IDX <- c()
for (x in rownames(g)) {

    idx <- which(df$phage_vOTU == x)

    if (length(idx) == 0) { next }

    IDX <- c(IDX, idx)
}
df <- df[IDX, ]



### write data
write.table(
    df,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'crAss_host_pairs.final.txt'
)
