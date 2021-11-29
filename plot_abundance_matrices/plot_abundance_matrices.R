.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(latex2exp)
sessionInfo()
source('function_plot_crass_tree.R')
source('function_plot_genera.R')
source('function_plot_sample_abundance.R')



######### READ DATA
tr <- read.tree('../from_Peregrine/MSA_crAss_TerL_MidpointRooted.2.newick')


DF <- read.table(
    '../from_Peregrine/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$cl_repres



# -------------------- data LLD --------------------
M1 <- read.table(
    '../from_Peregrine/LLD_crAss_abundance.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M1 <- M1[tr$tip.label, ]

k1 <- read.table(
    '../../LLD_LLD2_Info/LLD_GTMicrob.txt',
    sep = '\t',
    row.names = 1,
    header = FALSE,
    stringsAsFactors = FALSE
)

d1 <- read.table(
    '../../LLD_LLD2_Info/Pheno_science_imputed_1135/20150715_intristic_1135patients_log_imputed.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)
d1$work_id <- k1[rownames(d1), 1]
d1$work_id <- paste0('X', d1$work_id)



# -------------------- data LLD2 --------------------
M2 <- read.table(
    '../from_Peregrine/LLD2_crAss_abundance.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M2 <- M2[tr$tip.label, ]

k2 <- read.table(
    '../../LLD_LLD2_Info/key_LLD_baseline_fup_338sample_pairs.txt',
    sep = '\t',
    row.names = 3,
    header = TRUE,
    stringsAsFactors = FALSE
)

d2 <- read.table(
    '../../LLD_LLD2_Info/LLD2_Phenotypes/data_pheno_LLD_base_fup_338pairs_62pheno_18med_17disease_log_min_10participants.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)
d2$work_id <- k2[rownames(d2), 'barcode_fup']
d2$work_id <- sub('\\.bam$', '', d2$work_id)



# -------------------- data 300OB --------------------
M3 <- read.table(
    '../from_Peregrine/300OB_crAss_abundance.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M3 <- M3[tr$tip.label, ]

k3 <- read.table(
    '../../300OB_Info/key_300OB.txt',
    sep = '\t',
    header = TRUE,
    row.names = 6,
    stringsAsFactors = FALSE
)

d3 <- read.table(
    '../../300OB_Info/300OB_65phenotype.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
)
d3$work_id <- k3[rownames(d3), 'G_id']



# -------------------- data IBD --------------------
M4 <- read.table(
    '../from_Peregrine/IBD_crAss_abundance.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M4 <- M4[tr$tip.label, ]

d4 <- read.table(
    '../../IBD_Info/minimal_metadata_IBD_assembly.2.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)



######### CREATE PALETTE
V <- c(unlist(M1), unlist(M2), unlist(M3), unlist(M4))
V <- V[V != 0]
range <- range(V)
range <- log10(range)

br_power <- seq(from = floor(range[1]), to = ceiling(range[2]), by = 0.1)
brackets <- 10^br_power

colfunc <- colorRampPalette(c('skyblue1', 'skyblue2', 'red2', 'red1'))

P <- data.frame(
    col = colfunc(length(brackets) - 1),
    Above = head(brackets, -1),
    BelowEqual = tail(brackets, -1),
    stringsAsFactors = FALSE
)



######### PLOT
pdf('crAss_abundance_matrices.pdf', height = 8.75, width = 7.5)

N <- ncol(M1) + ncol(M2) + ncol(M3) + ncol(M4) + 40*4 + 30
layout(
    matrix(1:5, ncol = 5),
    width = c(0.8, 2*((ncol(M1)+40)/N), 2*((ncol(M2)+40)/N), 2*((ncol(M3)+40)/N), 2*((ncol(M4)+70)/N))
)

par(mar = c(0, 0, 5, 0), cex = 0.5)



# -------------------- plot tree --------------------
plot_crass_tree(
    tr,
    DF,
    xMax = 6.05,
    bootstrap = FALSE,
    tip_labels = FALSE,
    cex_group_labels = 3,
    adj_group_labels = list(delta = c(0.75, -3))
)

plot_genera(
    tr,
    DF,
    X1 = 5.5,
    X2 = 6,
    all_labels = FALSE,
    sele_labels = readRDS('../genera_abundance_plot/genera_lab.RDS'),
    sele_labels_2shift = 'alpha17'
)

text(x = 5.75, y = Ntip(tr) + 5, labels = 'Genus-level\nclusters', adj = 0, cex = 1.25, srt = 90, xpd = TRUE)

YLIM = par()$usr[3:4]



# -------------------- plot legend --------------------
X <- seq(from = 0.5, to = 3.5, length = length(br_power))
rect(
    xleft = X[-length(X)],
    xright = X[-1],
    ybottom = Ntip(tr) + 11,
    ytop =  Ntip(tr) + 20,
    col = P$col,
    border = NA,
    xpd = TRUE
)

br_power_lab <- br_power[ c(1, ceiling(length(br_power) / 2), length(br_power)) ]
text(
    labels = as.expression(sapply(br_power_lab, function(x) bquote(10^.(x)))),
    x = X[which(br_power %in% br_power_lab)],
    y = Ntip(tr) + 5,
    xpd = TRUE
)

text(
    labels = 'Abundance',
    x = (0.5 + 3.5) / 2,
    y = Ntip(tr) + 25,
    cex = 1.5,
    xpd = TRUE
)



# -------------------- plot LLD --------------------
d1 <- d1[d1$work_id %in% colnames(M1), ]

s_order <- d1$work_id[ order(d1$antrop_age) ]

M1 <- M1[, s_order]

plot(
    NA,
    xlim = c(-20, ncol(M1) + 20),
    ylim = YLIM,
    xaxs = 'i',  yaxs = 'i',
    ann = FALSE, axes = FALSE
)

invisible(outer(1:nrow(M1), 1:ncol(M1), plot_sample_abundance, matrix = M1, pal = P))

rect(xleft = 0, xright = ncol(M1) + 1, ybottom = 0, ytop = Ntip(tr) + 1)

text(
    x = ncol(M1) / 2,
    y = Ntip(tr) + 22,
    labels = paste0('LLD\n(n=', ncol(M1), ')'),
    cex = 1.5,
    xpd = TRUE
)



# -------------------- plot LLD2 --------------------
d2 <- d2[d2$work_id %in% colnames(M2), ]

s_order <- d2$work_id[ order(d2$antrop_age) ]

M2 <- M2[, s_order]

plot(
    NA,
    xlim = c(-20, ncol(M2) + 20),
    ylim = YLIM,
    xaxs = 'i',  yaxs = 'i',
    ann = FALSE, axes = FALSE
)

invisible(outer(1:nrow(M2), 1:ncol(M2), plot_sample_abundance, matrix = M2, pal = P))

rect(xleft = 0, xright = ncol(M2) + 1, ybottom = 0, ytop = Ntip(tr) + 1)

text(
    x = ncol(M2) / 2,
    y = Ntip(tr) + 18.5,
    labels = paste0('LLD\nfollow-up\n(n=', ncol(M2), ')'),
    cex = 1.5,
    xpd = TRUE
)



# -------------------- plot 300OB --------------------
d3 <- d3[d3$work_id %in% colnames(M3), ]

s_order <- d3$work_id[ order(d3$BMI) ]

M3 <- M3[, s_order]


plot(
    NA,
    xlim = c(-20, ncol(M3) + 20),
    ylim = YLIM,
    xaxs = 'i',  yaxs = 'i',
    ann = FALSE, axes = FALSE
)

invisible(outer(1:nrow(M3), 1:ncol(M3), plot_sample_abundance, matrix = M3, pal = P))

rect(xleft = 0, xright = ncol(M3) + 1, ybottom = 0, ytop = Ntip(tr) + 1)

text(
    x = ncol(M3) / 2,
    y = Ntip(tr) + 22,
    labels = paste0('300OB\n(n=', ncol(M3), ')'),
    cex = 1.5,
    xpd = TRUE
)



# -------------------- plot IBD --------------------
d4 <- d4[d4$work_id %in% colnames(M4), ]

d4$label <- ''
d4$label[d4$ibd_Diagnosis == 'CD' & d4$ibd_CurrentStomaOrPouch == 'no' ] <- 'CD'
d4$label[d4$ibd_Diagnosis == 'CD' & d4$ibd_CurrentStomaOrPouch == 'yes'] <- 'CD, stoma'
d4$label[d4$ibd_Diagnosis == 'UC' & d4$ibd_CurrentStomaOrPouch == 'no' ] <- 'UC'
d4$label[d4$ibd_Diagnosis == 'UC' & d4$ibd_CurrentStomaOrPouch == 'yes'] <- 'UC, stoma'

categor <- c('CD', 'CD, stoma', 'UC', 'UC, stoma', '')
s_order <- c()
axis_at <- c(0)

for (x in categor) {

    s <- d4$work_id[ d4$label == x ]
    s_order <- c(s_order, s)
    axis_at <- c(axis_at, tail(axis_at, 1) + length(s))

}

axis_at <- head(axis_at, -1)
M4 <- M4[, s_order]


plot(
    NA,
    xlim = c(-20, ncol(M4) + 50),
    ylim = YLIM,
    xaxs = 'i',  yaxs = 'i',
    ann = FALSE, axes = FALSE
)

invisible(outer(1:nrow(M4), 1:ncol(M4), plot_sample_abundance, matrix = M4, pal = P))

rect(xleft = 0, xright = ncol(M4) + 1, ybottom = 0, ytop = Ntip(tr) + 1)

text(
    x = ncol(M4) / 2,
    y = Ntip(tr) + 25,
    labels = paste0('IBD (n=', ncol(M4), ')'),
    cex = 1.5,
    xpd = TRUE
)


axis(side = 3, pos = Ntip(tr) + 1, at = axis_at, labels = FALSE)

invisible(sapply(1:4, function (i) {

    text(
        x = sum(axis_at[i:(i+1)]) / 2,
        y = Ntip(tr) + 5,
        labels = sub(', stoma', '*', categor[i]),
        adj = 0,
        cex = 1.5,
        srt = 90,
        xpd = TRUE
    )

}))


dev.off()
