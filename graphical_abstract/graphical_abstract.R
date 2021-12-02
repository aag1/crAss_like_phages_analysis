.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(vioplot)
source('../plot_abundance_matrices/function_plot_sample_abundance.R')
sessionInfo()




# -------------------------------------- Functions -------------------------------------- #
plot_person <- function (X = 0, Y = 0, H = 1, W = 1, LWD = 1, COL = NA) {

    a <- pi * seq(from = -0.5, to = 0.5, by = 0.05)

    Xcoo <- c(0, 1, 1, 4, 1, 0, cos(a))
    Ycoo <- c(0, -4, 2, 0, 3, 3, 4 + sin(a))

    polygon(
        x = W * c(Xcoo, -1 * rev(Xcoo)) + X,
        y = H * c(Ycoo, rev(Ycoo)) + Y,
        lwd = LWD,
        col = COL
    )

}


plot_dna <- function (X = 0, Y = 0, H = 1, W = 1, LWD = 1) {

    v <- pi * seq(from = -1, to = 1, by = 0.1)

    lines(
        x = 1 * W * cos(v) + X,
        y = H * v + Y,
        lend = 'butt',
        lwd = LWD
    )

    lines(
        x = -1 * W * cos(v) + X,
        y = H * v + Y,
        lend = 'butt',
        lwd = LWD
    )

    for (z in v[c(2,5,8,11,14,17,20)]) {

        lines(
            x = W * c(-1, 1) * cos(z) + X,
            y = H * rep(z, 2) + Y,
            lend = 'butt',
            lwd = LWD
        )

    }

}


plot_crAss <- function (X = 0, Y = 0, H = 1, W = 1, LWD = 1, COL = NA) {

    Xcoo <- c(0, 0.25, 0.5, 0)
    Ycoo <- c(-0.5, -0.5, 0.5, 0.5)

    polygon(
        x = W * c(Xcoo, -1 * rev(Xcoo)) + X,
        y = H * c(Ycoo, rev(Ycoo)) + Y,
        lwd = LWD,
        col = COL
    )

    Xcoo <- c(0, 2, 2, 0)
    Ycoo <- c(0, 1, 3, 4)

    polygon(
        x = W * c(Xcoo, -1 * rev(Xcoo)) + X,
        y = H * c(Ycoo, rev(Ycoo)) + Y,
        lwd = LWD,
        col = COL
    )

}




# -------------------------------------- Initiate plot -------------------------------------- #
pdf('graphical_abstract.pdf')

layout(
    matrix(c(1,1,1,2,3,4,5,5,6), nrow = 3, ncol = 3, byrow = T),
    width = c(1, 3, 4)
)




# -------------------------------------- Outline -------------------------------------- #
par(mar = rep(1, 4))

plot(NA, xlim = c(0, 100), ylim = c(-4, 12), ann = F, axes = F)

plot_person(X = 8, Y = 6.25, W = 1.25, LWD = 2, COL = 'grey92')
text('1950 individuals\n(healthy, obese,\nIBD patients)', x = 8, y = -2, cex = 1.6)

arrows(x0 = 23, x1 = 33, y0 = 7, col = 'steelblue3', lwd = 3)

plot_dna(X = 48, Y = 7, H = 1.5, W = 2.5, LWD = 2)
text('2291 fecal\nmetagenomes', x = 48, y = -2, cex = 1.6)

arrows(x0 = 63, x1 = 73, y0 = 7, col = 'steelblue3', lwd = 3)

plot_crAss(X = 88, Y = 4, LWD = 2, W = 2.5, H = 1.5, COL = 'grey92')
text('1556 crAss-like phage\ngenomes\n(TerL, portal, MCP)', x = 88, y = -2, cex = 1.6, xpd = T)




# -------------------------------------- Variation in abundance -------------------------------------- #
M <- read.table(
    '../from_Peregrine/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

M <- M[grep('\\|g__[^\\|]+$', rownames(M)), ]
rownames(M) <- sub('^.+\\|g__(.+)$', '\\1', rownames(M))

sele_rows <- which(sapply(1:nrow(M), function (i) (sum(M[i, ] > 0) / ncol(M)) > 0.05))
sele_cols <- 50 * 1:44
M <- M[sele_rows, sele_cols]



DF <- read.table(
    '../from_Peregrine/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

tr <- read.tree('../from_Peregrine/MSA_crAss_TerL_MidpointRooted.2.newick')

sele_tips <- c()
for (x in rownames(M)) { sele_tips <- c(sele_tips, DF$cl_repres[DF$genus == x][1]) }
tr <- keep.tip(tr, sele_tips)
tr <- read.tree(text = write.tree(tr))

tr$tip.label <- sapply(tr$tip.label, function (x) DF$genus[DF$cl_repres == x])
M <- M[tr$tip.label, ]



V <- c(unlist(M))
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



par(mar = c(2, 1.5, 6, 0))

plot(tr, show.tip.label = F)

abline(h = par()$usr[4] + par()$mai[3] * (diff(par()$usr[3:4]) / par()$pin[2]), col = 'steelblue3', lwd = 3, xpd = T)

YLIM = par()$usr[3:4]



par(mar = c(2, 0, 6, 4))

plot(
    NA,
    xlim = c(0, ncol(M) + 1),
    ylim = YLIM,
    xaxs = 'i',  yaxs = 'i',
    ann = FALSE, axes = FALSE
)

invisible(outer(1:nrow(M), 1:ncol(M), plot_sample_abundance, matrix = M, pal = P))

rect(xleft = 0.5, xright = ncol(M) + 0.5, ybottom = par()$usr[3], ytop = par()$usr[4], xpd = T)

text('samples', x = ncol(M) / 2, y = Ntip(tr) + 3.5, cex = 1.6, xpd = T)

text('vOTUs', x = ncol(M) + 3.5, y = sum(YLIM) / 2, cex = 1.6, xpd = T, srt = -90)

mtext('Variation in abundance', side = 3, line = 2.75, cex = 1.35, col = 'steelblue3', font = 2)

abline(h = par()$usr[4] + par()$mai[3] * (diff(par()$usr[3:4]) / par()$pin[2]), col = 'steelblue3', lwd = 3, xpd = T)




# -------------------------------------- Widespread recombination -------------------------------------- #
load('../from_Peregrine/delta27_OLNE01000081_simplot_data.RData')
PALETTE <- readRDS('../from_Peregrine/delta27_vOTUs_palette.rds')

par(mar = c(4, 5, 5, 1))

plot(
    NA,
    xlim = c(0, 110000),
    ylim = c(0, 100),
    xlab = 'reference genome',
    ylab = '% identity',
    axes = F, ann = F, bty = 'n'
)

rect(xleft = 48724, xright = 72821, ybottom = par()$usr[3], ytop = par()$usr[4], col = 'grey95', border = NA)
text('transcription\nmodule', x = (48724 + 72821) / 2, y = 12, cex = 1.2)

invisible(
    lapply(c('NL_crAss001078', 'MT774406.1', 'OHAY01000065'), function (x) {

        lines(x = L[[x]]$at, y = L[[x]]$sim, col = PALETTE[x], lwd = 1.5)
    })
)

lines(x = c(0, 105000), y = rep(par()$usr[3], 2), xpd = T)
text('reference genome', x = 50000, y = -25, cex = 1.6, xpd = T)

axis(side = 2, at = c(0, 100), las = 1, cex.axis = 1.6)
text('% identity', x = -22000, y = 45, cex = 1.6, xpd = T, srt = 90)

mtext('Widespread recombination', side = 3, line = 1.75, cex = 1.35, col = 'steelblue3', font = 2, at = 50000)

abline(v = par()$usr[1] - par()$mai[2] * (diff(par()$usr[1:2]) / par()$pin[1]), col = 'steelblue3', lwd = 3, xpd = T)
abline(h = par()$usr[4] + par()$mai[3] * (diff(par()$usr[3:4]) / par()$pin[2]), col = 'steelblue3', lwd = 3, xpd = T)




# -------------------------------------- 4-year stability -------------------------------------- #
load('../crAss_stability/crass_violin_plot_data.RData')

par(mar = c(5, 5.5, 5, 3))

vioplot(
    V3, V1,
    col = 'grey92',
    names = NA,
    frame.plot = F,
    ylim = c(0, 1), yaxs = 'i',
    xaxt = 'n', yaxt = 'n', ann = F
)

text(c('intra-', 'inter-'), x = c(1, 2), y = -0.2, cex = 1.6, xpd = T)
text('individual', x = 1.5, y = -0.4, cex = 1.6, xpd = T)

axis(side = 2, at = c(0, 1), las = 1, cex.axis = 1.6)
text('Dissimilarity', x = 0, y = 0.5, cex = 1.6, xpd = T, srt = 90)

mtext('4-year stability', side = 3, line = 1.75, cex = 1.35, col = 'steelblue3', font = 2)

abline(h = par()$usr[4] + par()$mai[3] * (diff(par()$usr[3:4]) / par()$pin[2]), col = 'steelblue3', lwd = 3, xpd = T)




# -------------------------------------- Links to human phenotypes -------------------------------------- #
par(mar = c(5, 6, 5, 2))

bp <- barplot(
        c(93.22, 69.45),
        col = 'grey92',
        ylim = c(0, 100),
        yaxs = 'i', yaxt = 'n',
        las = 1
)

text(c('general\npopulation', 'IBD\npatients'), x = bp, y = -30, cex = 1.6, xpd = T)

axis(side = 2, at = c(0, 100), las = 1, cex.axis = 1.6)
text('% crAss-positive', x = -0.41, y = 45, cex = 1.6, xpd = T, srt = 90)

mtext('Links to human phenotypes', side = 3, line = 1.75, cex = 1.35, col = 'steelblue3', font = 2, at = 1.2)
abline(v = par()$usr[1] - par()$mai[2] * (diff(par()$usr[1:2]) / par()$pin[1]), col = 'steelblue3', lwd = 3, xpd = T)
abline(h = par()$usr[4] + par()$mai[3] * (diff(par()$usr[3:4]) / par()$pin[2]), col = 'steelblue3', lwd = 3, xpd = T)

dev.off()
