.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
sessionInfo()



# -------------------- read data --------------------
DF <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)



# -------------------- merge cohort tables --------------------
TAB <- data.frame(NULL)

for (x in c('LLD', 'LLD2', '300OB', 'IBD')) {

    file <- paste0(x, '_crAss_abundance.txt')

    tab <- read.table(
        file,
        sep = '\t',
        header = TRUE,
        row.names = 1,stringsAsFactors = FALSE
    )

    if (ncol(TAB) == 0) { TAB <- tab } else { TAB <- cbind(TAB, tab) }

}

cat('\n\nRownames in DF and TAB are identical:', identical(sort(rownames(DF)), sort(rownames(TAB))), '\n\n\n')

TAB <- TAB[rownames(DF), ]



# -------------------- exclude non-gut crAss-like phages --------------------
non_gut <- which(DF[rownames(TAB), 'group'] == '')

cat('No non-gut crAss-like phages were detected:', all(TAB[non_gut, ] == 0), '\n\n\n')

TAB <- TAB[-c(non_gut), ]



# -------------------- add full taxonomy to species --------------------
SUPER <- TAB
rownames(SUPER) <- sapply(rownames(SUPER), function (x) paste0('o__crAss|f__', DF[x, 'group'], '|g__', DF[x, 'genus'], '|s__', x))



# -------------------- add genera --------------------
GENERA <- data.frame(NULL)

for (x in unique(DF$genus)) {

    if (x == '') { next }

    idx <- which(DF$genus == x)

    ids <- rownames(DF)[idx]

    v <- apply(TAB[ids, ], 2, sum)
    g <- as.data.frame(t(v), stringsAsFactors = FALSE)
    rownames(g) <- paste0('o__crAss|f__', DF$group[idx][1], '|g__', x)

    GENERA <- rbind(GENERA, g)

}

SUPER <- rbind(GENERA, SUPER)



# -------------------- add families --------------------
FAMILIES <- data.frame(NULL)

for (x in unique(DF$group)) {

    if (x == '') { next }

    idx <- which(DF$group == x)

    ids <- rownames(DF)[idx]

    v <- apply(TAB[ids, ], 2, sum)
    f <- as.data.frame(t(v), stringsAsFactors = FALSE)
    rownames(f) <- paste0('o__crAss|f__', x)

    FAMILIES <- rbind(FAMILIES, f)

}

SUPER <- rbind(FAMILIES, SUPER)



# -------------------- add order --------------------
v <- apply(TAB, 2, sum)
ORDER <- as.data.frame(t(v), stringsAsFactors = FALSE)
rownames(ORDER) <- 'o__crAss'

SUPER <- rbind(ORDER, SUPER)



# -------------------- write table --------------------
write.table(
    SUPER,
    sep = '\t',
    quote = FALSE,
    file = 'LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt'
)
