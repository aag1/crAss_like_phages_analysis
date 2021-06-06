.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(phangorn)
sessionInfo()



# -------------------- read data --------------------
df <- read.table('/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt', sep = '\t', header = T, stringsAsFactors = F)


PALETTE <- readRDS('delta27_vOTUs_palette.rds')


tr <- read.tree('delta27_RNAP_VP02740_MSA.fasta.treefile')

tr <- midpoint(tr)
tr <- read.tree(text = write.tree(tr))

for (n in c(218:221, 131:132, 163, 117, 135)) {
    tr <- rotate(tr, node = n)
    tr <- read.tree(text = write.tree(tr))
}



# -------------------- PLOT --------------------
pdf('delta27_RNAP_VP02740_tree.pdf', height = 8.75, width = 3.5)

par(mar = c(0, 0.5, 0, 0.5), cex = 0.8)


# tree
plot.phylo(tr, x.lim = 4.6, show.tip.label = FALSE)

add.scale.bar()


# colored dots
genomes <- sub('_RNAP_VP02740$', '', tr$tip.label)

vOTUs   <- sapply(genomes, function (x) df$cl_repres[df$cl_member == x])

tiplabels(
    pch = 20,
    col = PALETTE[ vOTUs ],
    offset = 0.1
)


# tip labels
for (i in which(genomes == vOTUs)) {

    ADJ <- c(0, 0.5)
    if (genomes[i] == 'cs_ms_40_1') { ADJ <- c(0, 0) }
    if (genomes[i] == 'NL_crAss000239') { ADJ <- c(0, 1) }

    tiplabels(
        genomes[i],
        tip = i,
        bg = NULL,
        frame = 'none',
        adj = ADJ,
        offset = 0.25,
        xpd = TRUE
    )

}


# bootstrap
d <- as.data.frame(tr$edge, stringsAsFactors = FALSE)
colnames(d) <- c('nodeA', 'nodeD')
d$edge_len <- tr$edge.length
d <- d[ d$nodeD > Ntip(tr), ]
d <- d[ d$edge_len > 0.1, ]

nodelabels(
    tr$node.label[d$nodeD - Ntip(tr)],
    node = d$nodeD,
    bg = NULL,
    frame = 'none',
    adj = c(1.25, -0.25)
)


dev.off()
