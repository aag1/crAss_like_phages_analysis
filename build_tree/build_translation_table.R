sessionInfo()



### read data
df <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl_summary.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE)
rownames(df) <- df$cl_repres

tab <- read.table(
    '/data/umcg-tifn/DATABASES/data_Yutin_2020/TableS1_genomes_used_1.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE)

v1 <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/build_tree/TAGq.txt',
    header = FALSE,
    stringsAsFactors = FALSE)[, 1]

v2 <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/build_tree/TGAw.txt',
    header = FALSE,
    stringsAsFactors = FALSE)[, 1]



### add column
df$translation <- 'c11'

idx <- which(rownames(tab) %in% rownames(df))
df[ rownames(tab)[idx], 'translation' ] <- tab$translation[idx]
df[ 'OFRY01000050', 'translation' ] <- 'TAG|q'
df$translation <- sub('GenBank \\(c11\\)', 'c11', df$translation)

df[ v1, 'translation' ] <- 'TAG|q'

df[ v2, 'translation' ] <- 'TGA|w'



### write table
write.table(
    df,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'CRASS_DB_cl_summary_TRANSL.txt'
)
