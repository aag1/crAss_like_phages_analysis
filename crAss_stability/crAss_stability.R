library(latex2exp)
library(vegan)
library(vioplot)
sessionInfo()



# -------------------- read data --------------------
DF <- read.table(
    'CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$cl_repres


M <- read.table(
    'LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M <- M[grep('\\|g__[^\\|]+$', rownames(M)), ]
rownames(M) <- sub('^.+\\|g__(.+)$', '\\1', rownames(M))


t <- read.table(
    'sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)



# -------------------- data LLD --------------------
M1 <- M[, paste0('X', t$sample_id[t$cohort == 'LLD'])]

k1 <- read.table(
    'LLD_GTMicrob.txt',
    sep = '\t',
    row.names = 1,
    header = FALSE,
    stringsAsFactors = FALSE
)

d1 <- read.table(
    '20150715_intristic_1135patients_log_imputed.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)
d1$work_id <- k1[rownames(d1), 1]
d1$work_id <- paste0('X', d1$work_id)



# -------------------- data LLD2 --------------------
M2 <- M[, t$sample_id[t$cohort == 'LLD2']]

k2 <- read.table(
    'key_LLD_baseline_fup_338sample_pairs.txt',
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

d1 <- d1[d1$work_id %in% colnames(M2), ]

s_order <- d1$work_id[ order(d1$antrop_age) ]

M1 <- M1[, s_order]
M2 <- M2[, s_order]



# -------------------- plot 1: stability of the individual phages --------------------
DF1 <- data.frame(
    LLD     = sapply(rownames(M1), function (x) sum(M1[x,] > 0)),
    LLD2    = sapply(rownames(M1), function (x) sum(M2[x,] > 0)),
    total   = sapply(rownames(M1), function (x) sum((M1[x,] > 0) | (M2[x,] > 0))),
    lost    = sapply(rownames(M1), function (x) sum((M1[x,] > 0) & (M2[x,] == 0))),
    aquired = sapply(rownames(M1), function (x) sum((M1[x,] == 0) & (M2[x,] > 0))),
    stable  = sapply(rownames(M1), function (x) sum((M1[x,] > 0) & (M2[x,] > 0))),
    stringsAsFactors = FALSE
)


sele <- which((DF1$LLD > ncol(M1)*0.1) | (DF1$LLD2 > ncol(M1)*0.1))
DF1 <- DF1[sele, ]  # select phages present in >10% samples


for (x in c('lost', 'aquired', 'stable')) { DF1[, x] <- DF1[, x] / DF1$total * 100 }
df1 <- DF1[, c('lost', 'aquired', 'stable')]


pdf('crAss_stability_1.pdf', height = 3.3, width = 7.5)
layout(matrix(1:2, ncol = 2))


Yaxis2 <- pretty(c(0, max(DF1$total)))
par(mar = c(2.25, 3.5, 2.75, 3.5), mgp = c(3, 0.7, 0), ps = 9)


bp <- barplot(
    t(as.matrix(df1)),
    space = 0,
    border = NA,
    col = c('forestgreen', 'lightcoral', 'grey90'),
    names = rep('', nrow(df1)),
    axes = FALSE, ann = FALSE
)

for (i in 1:nrow(df1)) {

    x <- rownames(df1)[i]
    g <- sub('[0-9]+$', '', x)
    n <- sub('^[a-z]+', '', x)
    text(x = bp[i], y = -5, labels = TeX(paste0('$\\', g, '$', n)), adj = 1, srt = 90, xpd = TRUE)

}

lines(x = bp, y = DF1$total / max(Yaxis2) * 100, lwd = 1.5)

axis(side = 2, pos = bp[1] - 0.5, las = 1)
axis(side = 4, pos = bp[length(bp)] + 0.5, at = Yaxis2 / max(Yaxis2) * 100, labels = Yaxis2, las = 1)

mtext(side = 2, line = 1.75, text = 'Total positive individuals, %', cex = 10/9)
text(x = par()$usr[2] + diff(par()$usr[1:2])*0.2, y = 50, labels = 'Total positive individuals', srt = -90, xpd = TRUE, cex = 10/9)

legend(
    x = par()$usr[1] + diff(par()$usr[1:2])*0.5,
    y = 118,
    fill = c('forestgreen', 'lightcoral', 'grey90'),
    legend = c('Baseline', 'Follow-up', 'Both'),
    ncol = 3,
    cex = 0.72,
    xjust = 0.5, xpd = TRUE,
)

mtext('A', side = 3, line = 1, at = -3.8, family = 'sans', cex = 12/9 * 2, xpd = TRUE)



# -------------------- numbers: stability of the individual phages --------------------
mean_stable <- mean(df1$stable)
cat('\nMean % stable samples:', mean_stable, '\n')

idx <- which(df1$stable > 50)
cat('Phages stable in >50% samples:', paste(rownames(df1)[idx], collapse = ', '), '\n\n')



# -------------------- plot 1: stability of the overal phage composition --------------------
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


par(mar = c(2.25, 6.5, 2.75, 0.5))

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

mtext('B', side = 3, line = 1, at = -0.5, family = 'sans', cex = 12/9 * 2, xpd = TRUE)

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



# -------------------- plot 2: stability in individuals --------------------
DF2 <- data.frame(
    total   = sapply(colnames(M1), function (x) sum((M1[,x] > 0) | (M2[,x] > 0))),
    lost    = sapply(colnames(M1), function (x) sum((M1[,x] > 0) & (M2[,x] == 0))),
    aquired = sapply(colnames(M1), function (x) sum((M1[,x] == 0) & (M2[,x] > 0))),
    stable  = sapply(colnames(M1), function (x) sum((M1[,x] > 0) & (M2[,x] > 0))),
    stringsAsFactors = FALSE
)


sele <- which(DF2$total > 0)
DF2 <- DF2[sele, ]
DF2 <- DF2[order(DF2$total, DF2$stable, DF2$aquired, decreasing = TRUE), ]


for (x in c('lost', 'aquired', 'stable')) { DF2[, x] <- DF2[, x] / DF2$total * 100 }
df2 <- DF2[, c('lost', 'aquired', 'stable')]


pdf('crAss_stability_2.pdf', width = 3.3, height = 3.3)


Yaxis2 <- pretty(c(0, max(DF2$total)))
par(mar = c(2.25, 3, 2.75, 2.5), mgp = c(3, 0.7, 0), ps = 9)


bp <- barplot(
    t(as.matrix(df2)),
    space = 0,
    border = NA,
    col = c('forestgreen', 'lightcoral', 'grey90'),
    names = rep('', nrow(df2)),
    axes = FALSE, ann = FALSE
)

lines(x = bp, y = DF2$total / max(Yaxis2) * 100)

axis(side = 2, pos = bp[1] - 0.5, las = 1)
axis(side = 4, pos = bp[length(bp)] + 0.5, at = Yaxis2 / max(Yaxis2) * 100, labels = Yaxis2, las = 1)

mtext(side = 1, line = 0.75, text = paste(length(sele), 'individuals'), cex = 10/9)
mtext(side = 2, line = 1.75, text = 'Total phages detected, %', cex = 10/9)
text(x = par()$usr[2] + diff(par()$usr[1:2])*0.15, y = 50, labels = 'Total phages detected', srt = -90, xpd = TRUE, cex = 10/9)

legend(
    x = par()$usr[1] + diff(par()$usr[1:2])*0.5,
    y = 118,
    fill = c('forestgreen', 'lightcoral', 'grey90'),
    legend = c('Baseline', 'Follow-up', 'Both'),
    ncol = 3,
    cex = 0.72,
    xjust = 0.5, xpd = TRUE,
)

dev.off()
