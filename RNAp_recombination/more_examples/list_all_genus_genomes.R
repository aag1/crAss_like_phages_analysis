sessionInfo()


genus <- commandArgs(trailingOnly = TRUE)[1]


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


genus_cl_repres <- t2$cl_repres[t2$genus == genus]

genus_all_genomes <- t1$cl_member[t1$cl_repres %in% genus_cl_repres]


write.table(
    genus_all_genomes,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    file = paste0(genus, '_genomes_all.ids')
)
