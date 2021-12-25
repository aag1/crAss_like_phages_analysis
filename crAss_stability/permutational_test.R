.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(vegan)
sessionInfo()



# -------------------- read data --------------------
DF <- read.table(
    '../from_Peregrine/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$cl_repres


M <- read.table(
    '../from_Peregrine/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M <- M[grep('\\|g__[^\\|]+$', rownames(M)), ]
rownames(M) <- sub('^.+\\|g__(.+)$', '\\1', rownames(M))


t <- read.table(
    '../from_Peregrine/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)



# -------------------- data LLD --------------------
M1 <- M[, paste0('X', t$sample_id[t$cohort == 'LLD'])]

k1 <- read.table(
    '../../LLD_LLD2_Info/LLD_GTMicrob.txt',
    sep = '\t',
    row.names = 1,
    header = FALSE,
    stringsAsFactors = FALSE
)



# -------------------- data LLD2 --------------------
M2 <- M[, t$sample_id[t$cohort == 'LLD2']]

k2 <- read.table(
    '../../LLD_LLD2_Info/key_LLD_baseline_fup_338sample_pairs.txt',
    sep = '\t',
    row.names = 3,
    header = TRUE,
    stringsAsFactors = FALSE
)



# -------------------- match LLD & LLD2 data --------------------
MATCH <- sapply(colnames(M2), function (lld2_work_id) {

    lld_work_id <- k2$barcode_baseline[ k2$barcode_fup == paste0(lld2_work_id, '.bam') ]

    lld_work_id <- sub('^fece_', '', lld_work_id)

    lld_work_id <- paste0('X', lld_work_id)

    if (lld_work_id == 'ZZZZ') { lld_work_id <- 'QQQQ' }

    return(lld_work_id)

})

colnames(M2) <- MATCH
M1 <- M1[, MATCH]



# -------------------- test stability of the overal phage composition --------------------
# Wilcoxon test on Bray-Curtis dissimilarity, with empiric P-value based on 10,000 permutations
# See Chen et al. 2021 (https://doi.org/10.1016/j.cell.2021.03.024)

colnames(M2) <- paste0(colnames(M2), '_F')
m <- cbind(M1, M2)

m <- t(m)
m <- m[, colSums(m)>0]
m <- m[rowSums(m)>0, ]

bl_ids <- grep('_F', row.names(m), value = TRUE, invert = TRUE)    # baseline sample ids
fu_ids <- grep('_F', row.names(m), value = TRUE, invert = FALSE)   # follow-up sample ids


# 10,000 permutations
PERM_RES <- data.frame(NULL)

for (i in 1:10000) {

    p <- m
    row.names(p) <- sample(row.names(m), nrow(m))


    D <- vegdist(p, method = 'bray')
    D <- as.matrix(D)


    d1 <- D[bl_ids, bl_ids]
    V1 <- d1[upper.tri(d1)]

    d2 <- D[fu_ids, fu_ids]
    V2 <- d2[upper.tri(d2)]

    sele <- which(paste0(bl_ids, '_F') %in% fu_ids)
    d3 <- D[bl_ids[sele], paste0(bl_ids[sele], '_F')]
    V3 <- diag(d3)


    res <- data.frame(
        pval_bl_fu = wilcox.test(V1, V2, alternative = 'two.sided', paired = F)$p.value,
        pval_in_bl = wilcox.test(V3, V1, alternative = 'two.sided', paired = F)$p.value,
        pval_in_fu = wilcox.test(V3, V2, alternative = 'two.sided', paired = F)$p.value,
        stringsAsFactors = FALSE
    )
    PERM_RES <- rbind(PERM_RES, res)


    if (i %% 100 == 0) { cat('Iteration', i, '...\n') }

}


# real data
D <- vegdist(m, method = 'bray')
D <- as.matrix(D)


d1 <- D[bl_ids, bl_ids]
V1 <- d1[upper.tri(d1)]

d2 <- D[fu_ids, fu_ids]
V2 <- d2[upper.tri(d2)]

sele <- which(paste0(bl_ids, '_F') %in% fu_ids)
d3 <- D[bl_ids[sele], paste0(bl_ids[sele], '_F')]
V3 <- diag(d3)


REAL_RES <- data.frame(
    pval_bl_fu = wilcox.test(V1, V2, alternative = 'two.sided', paired = F)$p.value,
    pval_in_bl = wilcox.test(V3, V1, alternative = 'two.sided', paired = F)$p.value,
    pval_in_fu = wilcox.test(V3, V2, alternative = 'two.sided', paired = F)$p.value,
    stringsAsFactors = FALSE
)


# empiric P-value
emp_pval_bl_fu <- sum( PERM_RES$pval_bl_fu < REAL_RES$pval_bl_fu ) / 10000
emp_pval_in_bl <- sum( PERM_RES$pval_in_bl < REAL_RES$pval_in_bl ) / 10000
emp_pval_in_fu <- sum( PERM_RES$pval_in_fu < REAL_RES$pval_in_fu ) / 10000


# save test results
save(PERM_RES, REAL_RES, emp_pval_bl_fu, emp_pval_in_bl, emp_pval_in_fu, file = 'crass_permutation_test_data.RData')
