.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(latex2exp)
sessionInfo()
source('/home/umcg-agulyaeva/crAss_analysis/build_tree/function_plot_crass_tree.R')
source('/home/umcg-agulyaeva/crAss_analysis/build_tree/function_plot_genera.R')



######### READ DATA
tr <- read.tree('/data/umcg-tifn/crAss_analysis/build_tree/MSA_crAss_TerL_MidpointRooted.2.newick')


DF <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$cl_repres


tab <- read.table(
    'crAss_host_pairs.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
tab <- tab[tab$phage_vOTU %in% tr$tip.label, ]



######### PLOT
pdf('crAss_host_crispr.pdf', height = 8.75, width = 7.5/2.8*2.3)

layout(
    matrix(1:3, ncol = 3),
    width = c(0.8, 0.5, 1)
)

par(mar = c(0, 0, 5, 0), cex = 0.5)



# -------------------- plot crAss tree --------------------
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
    sele_labels = unique(tab$phage_genus),
    sele_labels_2shift = c('alpha20', 'delta34', 'zeta5', 'zeta11', 'zeta24')
)

text(x = 5.75, y = Ntip(tr) + 5, labels = 'Genus-level\nclusters', adj = 0, cex = 1, srt = 90, xpd = TRUE)

text('crAss-like\nphages', x = sum(par()$usr[1:2]) / 2, y = par()$usr[4] + diff(par()$usr[3:4]) * 0.02, cex = 2, xpd = TRUE)

YLIM = par()$usr[3:4]



# -------------------- plot connections--------------------
tab$host_genus <- unlist(lapply(strsplit(tab$host_taxonomy, ';'), function (v) {

    x <- sub('^uncultured ', '', rev(v)[1])
    y <- strsplit(x, ' ')[[1]][1]
    return(y)

}))

tab <- tab[c(which(tab$phage_group == 'epsilon'), which(tab$phage_group != 'epsilon')), ]

tab <- tab[, c('phage_vOTU', 'host_genus')]

tab <- unique(tab)


V0 <- unique(tab$host_genus)
V <- rev(c('Fusobacterium', 'Runella', 'Odoribacter', 'Bacteroides', 'Prevotella', 'Porphyromonas', 'Parabacteroides', 'Bacillus', 'Streptococcus', 'Clostridium'))
if (!all(V0 %in% V)) { stop('Some host genera are missing from the plot: ', paste(V0[!(V0 %in% V)], collapse = ','), '!') }

K <- diff(YLIM) / length(V)


PALETTE <- setNames(
    c('#AA4488', '#4477AA', '#44AA77', '#AAAA44', '#AA7744'),
    c('alpha_gamma', 'beta', 'delta', 'epsilon', 'zeta')
)


plot(NA, xlim = c(0, 1), ylim = YLIM, xaxs = 'i', yaxs = 'i', ann = FALSE, axes = FALSE)

text('CRISPR\nlink', x = 0.5, y = par()$usr[4] + diff(par()$usr[3:4]) * 0.02, cex = 2, xpd = TRUE)

for (i in 1:nrow(tab)) {

    Y1 <- which(tr$tip.label == tab$phage_vOTU[i])

    Y2 <- YLIM[1] + (which(V == tab$host_genus[i]) - 0.5) * K

    group <- DF[tab$phage_vOTU[i], 'group']

    lines(
        x = c(0, 1),
        y = c(Y1, Y2),
        col = ifelse(group == '', 'black', PALETTE[ group ])
    )

}



# -------------------- plot hosts --------------------
plot(NA, xlim = c(0, 1), ylim = c(0.5, length(V) + 0.5), xaxs = 'i', yaxs = 'i', ann = FALSE, axes = FALSE)

text('Predicted\nhosts', x = 0.5, y = par()$usr[4] + diff(par()$usr[3:4]) * 0.02, cex = 2, xpd = TRUE)

text(
    x = 0.1,
    y = seq_along(V),
    labels = V,
    font = 3,
    adj = 0,
    cex = 2
)

lines(x = rep(0.65, 2), y = c(4, 8) + c(-0.25, 0.25))
text('Bacteroidales', x = 0.7, y = 6, srt = 90, cex = 2)

lines(x = rep(0.85, 2), y = c(4, 9) + c(-0.25, 0.25))
text('Bacteroidetes', x = 0.9, y = 6.5, srt = 90, cex = 2)

lines(x = rep(0.85, 2), y = c(1, 3) + c(-0.25, 0.25))
text('Firmicutes', x = 0.9, y = 2, srt = 90, cex = 2)



dev.off()
