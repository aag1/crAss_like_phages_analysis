sessionInfo()


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


delta27_cl_repres <- t2$cl_repres[t2$genus == 'delta27']

delta27_all_genomes <- t1$cl_member[t1$cl_repres %in% delta27_cl_repres]


write.table(
    delta27_all_genomes,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = 'delta27_genomes_all.ids'
)
