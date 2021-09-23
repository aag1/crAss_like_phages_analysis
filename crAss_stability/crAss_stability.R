.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(latex2exp)
library(vegan)
library(vioplot)
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



# -------------------- plot stability of the individual phages --------------------
DF1 <- data.frame(
    LLD_pos  = sapply(rownames(M1), function (x) sum(M1[x,] > 0)),
    LLD2_pos = sapply(rownames(M1), function (x) sum(M2[x,] > 0)),
    both_pos_obs = sapply(rownames(M1), function (x) sum((M1[x,] > 0) & (M2[x,] > 0))),
    stringsAsFactors = FALSE
)
DF1$both_pos_exp <- ncol(M1) * (DF1$LLD_pos / ncol(M1)) * (DF1$LLD2_pos / ncol(M1))


sele <- which((DF1$LLD_pos > ncol(M1)*0.1) | (DF1$LLD2_pos > ncol(M1)*0.1))
DF1 <- DF1[sele, ]  # select phages present in >10% samples


pdf('crAss_stability.pdf', height = 3.3, width = 7.5)
layout(matrix(1:2, ncol = 2))

par(mar = c(2.25, 4.5, 2.75, 2), ps = 9)

plot(
    NA,
    xlim = c(0, nrow(DF1)),
    ylim = c(0, 100),
    xaxs = 'i', yaxs = 'i',
    xaxt = 'n',
    las = 1,
    bty = 'n',
    ann = FALSE
)

mtext(side = 2, line = 2.25, text = 'Individuals positive in both timepoints', cex = 10/9)

points(
    x = 1:nrow(DF1),
    y = DF1$both_pos_exp,
    pch = 4,
    cex = 0.75,
    xpd = TRUE
)

points(
    x = 1:nrow(DF1),
    y = DF1$both_pos_obs,
    pch = 16,
    cex = 0.75,
    xpd = TRUE
)

for (i in 1:nrow(DF1)) {

    x <- rownames(DF1)[i]
    g <- sub('[0-9]+$', '', x)
    n <- sub('^[a-z]+', '', x)
    text(x = i, y = -5, labels = TeX(paste0('$\\', g, '$', n)), adj = 1, srt = 90, xpd = TRUE)

}

legend(
    'top',
    inset = -0.2,
    pch = c(4, 16),
    pt.cex = 0.75,
    legend = c('expected', 'observed'),
    ncol = 2,
    xpd = TRUE
)

mtext('A', side = 3, line = 1, at = -3.8, family = 'sans', cex = 12/9 * 2, xpd = TRUE)



# -------------------- plot stability of the overal phage composition --------------------
m <- t(M)
m <- m[,colSums(m)>0]
m <- m[rowSums(m)>0,]

D <- vegdist(m, method = 'bray')
D <- as.matrix(D)


sele <- MATCH
sele <- sele[sele %in% rownames(D)]
d1 <- D[sele, sele]
V1 <- d1[lower.tri(d1)]

sele <- names(MATCH)
sele <- sele[sele %in% rownames(D)]
d2 <- D[sele, sele]
V2 <- d2[lower.tri(d2)]

MATCH <- MATCH[MATCH %in% rownames(D)]
MATCH <- MATCH[names(MATCH) %in% rownames(D)]
V3 <- sapply(seq_along(MATCH), function (i) D[MATCH[i], names(MATCH)[i]])

cat('\n# datapoints for violinl plot:', length(V3), length(V1), length(V2), '\n\n')


W1 <- wilcox.test(V3, V1, alternative = 'less', paired = FALSE)
W2 <- wilcox.test(V3, V2, alternative = 'less', paired = FALSE)


par(mar = c(2.25, 5.5, 2.75, 0.5))

vioplot(
    V3, V1, V2,
    col = c('grey60', 'forestgreen', 'lightcoral'),
    names = NA,
    las = 1,
    frame.plot = FALSE,
    xaxt = 'n', bty = 'n'
)

text(
    x = c(1, 2.15, 3.15),
    y = -0.12,
    labels = c('Inter-individual', 'Baseline', 'Follow-up'),
    xpd = TRUE
)

mtext(side = 2, line = 2.25, text = 'Bray-Curtis distance', cex = 10/9)

lines(x = c(1, 2), y = rep(1.06, 2), xpd = TRUE)
text(x = 1.5, y = 1.09, labels = formatC(W1$p.value, format = 'e', digits = 2), cex = 8/9, xpd = TRUE)

lines(x = c(1, 3), y = rep(1.18, 2), xpd = TRUE)
text(x = 2, y = 1.21, labels = formatC(W2$p.value, format = 'e', digits = 2), cex = 8/9, xpd = TRUE)

mtext('B', side = 3, line = 1, at = -0.55, family = 'sans', cex = 12/9 * 2, xpd = TRUE)

dev.off()


cat('\nMean intra-individual BC:', mean(V3), '\n')
cat('\nMean inter-individual BC LLD:', mean(V1), '\n')
cat('\nMean inter-individual BC FUP:', mean(V2), '\n\n')
cat('\nIntra-individual BC vs. LLD BC p-value:', W1$p.value, '\n')
cat('\nIntra-individual BC vs. FUP BC p-value:', W2$p.value, '\n\n')


# inquiry into low inter-individual BC distances
pdf('boxplots_test.pdf')

boxplot(
    V3, V1, V2,
    col = c('grey60', 'forestgreen', 'lightcoral'),
    names = NA,
    las = 1,
    frame.plot = FALSE,
    xaxt = 'n', bty = 'n'
)

V3ext <- c(V3, rep(1, 338 - length(V3)))
boxplot(
    V3ext, V1, V2,
    col = c('grey60', 'forestgreen', 'lightcoral'),
    names = NA,
    las = 1,
    frame.plot = FALSE,
    xaxt = 'n', bty = 'n'
)

dev.off()


cat ('\nMost close LLD samples:\n')

min(V1)

samples <- rownames(which(d1 == min(V1), arr.ind = TRUE))

taxa <- rownames(M)[ (M[, samples[1]] > 0) | (M[, samples[2]] > 0) ]

M[taxa, samples, drop = FALSE]


cat ('\nMost close FUP samples:\n')

min(V2)

samples <- rownames(which(d2 == min(V2), arr.ind = TRUE))

taxa <- rownames(M)[ (M[, samples[1]] > 0) | (M[, samples[2]] > 0) ]

M[taxa, samples, drop = FALSE]

cat('\n\n')
