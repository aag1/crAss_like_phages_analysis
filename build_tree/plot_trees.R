.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(phangorn)
library(latex2exp)
sessionInfo()
source('/home/umcg-agulyaeva/crAss_analysis/build_tree/function_plot_crass_tree.R')
source('/home/umcg-agulyaeva/crAss_analysis/build_tree/function_plot_tree_support_pch.R')
source('/home/umcg-agulyaeva/crAss_analysis/build_tree/function_plot_genera.R')



### read data
DF <- read.table(
    'CRASS_DB_cl_summary_TRANSL_GROUPS.1.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$cl_repres



################################################## TerL tree ##################################################
tr <- read.tree('MSA_crAss_TerL_MidpointRooted.1.newick')



### rotate tree
for (n in c(372, 495, 556, 670, 682, 684, 689, 698, 637, 638, 651:653, 557, 558, 496, 539, 373:374, 395, 399, 397, 396, 407:436, 449:452, 478, 656, 437, 438, 443, 444)) {
    tr <- rotate(tr, node = n)
    tr <- read.tree(text = write.tree(tr))
}

write.tree(tr, file = 'MSA_crAss_TerL_MidpointRooted.2.newick')



### order table by the tree, introduce nice genera names
non_tree <- rownames(DF)[!(rownames(DF) %in% tr$tip.label)]

DF <- DF[c(rev(tr$tip.label), non_tree), ]


V <- setNames(rep(0, 6), c('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta'))

GAMMA <- c('NL_crAss001398', 'OATI01000137', 'OLND01000117', 'OGRD01000043', 'LSPZ01000006', 'HvCF_A12_ms_2', 'NL_crAss000735', 'NL_crAss000806', 'OLGH01000138')

DF$genus <- ''

for (x in unique(DF$genus_raw)) {

    idx <- which(DF$genus_raw == x)

    group <- unique(DF$group[idx])


    if (length(group) > 1) { stop('Members of multiple groups in one genus: ', x, '!') }

    if (group == 'alpha_gamma') { group <- ifelse(all(rownames(DF)[idx] %in% GAMMA), 'gamma', 'alpha') }

    if (!(group %in% names(V))) { next }


    N <- V[group] + 1

    genus <- paste0(group, N)

    DF$genus[idx] <- genus

    V[group] <- N

}

DF$genus_raw <- NULL


write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'CRASS_DB_cl_summary_TRANSL_GROUPS.2.txt'
)


DF$n_this_study <- DF$n_LLD + DF$n_LLD2 + DF$n_300OB + DF$n_IBD

DF$n_DBs <- DF$n_db673 + DF$n_db249 + DF$n_db146

DF <- DF[, c('cl_repres', 'cl_size', 'n_this_study', 'n_DBs', 'translation', 'group', 'genus')]

write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt'
)



### plot tree
pdf('crAss_tree_TerL.pdf', height = 11.7*2, width = 8.3)

par(mar = rep(0.5, 4))

plot_crass_tree(
    tr,
    DF,
    xMax = 5.5,
    bootstrap = TRUE,
    tip_labels = TRUE,
    cex_tip_labels = 0.3,
    cex_group_labels = 3
)

#nodelabels(col = 'red', frame = 'none', adj = c(1.1, -0.5))

plot_genera(tr, DF, X1 = 5, X2 = 5.5, all_labels = TRUE, cex_labels = 0.3)

text(x = 5.25, y = Ntip(tr) + 5, labels = 'Genus-level\nclusters')

dev.off()



################################################## portal tree ##################################################
tr <- read.tree('MSA_crAss_portal.fasta.treefile')



### root tree
tr <- midpoint(tr)
tr <- read.tree(text = write.tree(tr))



### rotate tree
for (n in c(385, 677, 686, 688, 689, 565, 448, 449, 394)) {
    tr <- rotate(tr, node = n)
    tr <- read.tree(text = write.tree(tr))
}

write.tree(tr, file = 'MSA_crAss_portal_MidpointRooted.newick')



### plot tree
pdf('crAss_tree_portal.pdf', height = 11.7*2, width = 8.3)

par(mar = rep(0.5, 4))

plot_crass_tree(
    tr,
    DF,
    xMax = 6.25,
    bootstrap = TRUE,
    tip_labels = TRUE,
    cex_tip_labels = 0.3,
    cex_group_labels = 3,
    adj_group_labels = list(beta = c(1.5, 1.5))
)

#nodelabels(col = 'red', frame = 'none', adj = c(1.1, -0.5))

plot_genera(tr, DF, X1 = 5.75, X2 = 6.25, all_labels = TRUE, cex_labels = 0.3)

text(x = 6, y = Ntip(tr) + 5, labels = 'Genus-level\nclusters')

dev.off()
