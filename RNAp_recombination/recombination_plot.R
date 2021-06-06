.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(bio3d)
source('/home/umcg-agulyaeva/crAss_analysis/genome_maps/function_plot_genome_map.R')
sessionInfo()




# -------------------- read data --------------------
R <- readRDS('/data/umcg-tifn/crAss_analysis/map_reads/LLD_crAss_depth.rds')


d <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/genome_maps/selected_domains.txt',
    sep = '\t',
    header = FALSE,
    comment.char = '',
    stringsAsFactors = FALSE)

COL <- setNames(d[, 2], d[, 1])
COL <- c(COL, VP02740 = 'navajowhite1')


tabF <- '/data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/OLNE01000081_TAGq_AA_coordinates.txt'
tab <- read.table(tabF, sep = '\t', header = FALSE, stringsAsFactors = FALSE)


domF <- '/data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/crAss_key_doms_in_OLNE01000081_TAGq.txt'
dom <- read.table(domF, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

dom <- dom[rownames(dom) != 'OLNE01000081_12',]
dom['OLNE01000081_11', 'TerL'] <- '1-833'
dom['OLNE01000081_39', 'RNApBprime'] <- '1-5441'
dom['OLNE01000081_40',] <- '-'
dom[, 'VP02740'] <- '-'
dom['OLNE01000081_40', 'VP02740'] <- '1-1945'   # VP02740 nickname is 'long_transcr'


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


MSA <- read.fasta('delta27_genomes_circular_MSA.fasta')$ali




# -------------------- calculations for simplot --------------------
sele_seq <- 'OLNE01000081'
other_seq <- rownames(MSA)[rownames(MSA) != sele_seq]


idx <- which(MSA[sele_seq, ] != '-')
MSA <- MSA[, idx]


w_center <- seq(from = 501, to = ncol(MSA) - 500, by = 200)
w_from <- w_center - 500
w_to <- w_center + 500
w_to[ length(w_to) ] <- ncol(MSA)
w_len <- w_to - w_from + 1


L <- lapply(other_seq, function (X) {

    ali <- MSA[c(X, sele_seq), ]


    df <- data.frame(
        from = w_from,
        to   = w_to,
        len  = w_len,
        at   = w_center,
        stringsAsFactors = FALSE
    )


    # Four genomes that are close to completeness and represent vOTUs were included as well.
    if (X %in% c('cs_ms_40_1', 'NL_crAss001084', 'NL_crAss001078', 'NL_crAss000719')) {

        coo <- which(ali[X,] != '-')
        df <- df[df$from >= min(coo), ]
        df <- df[df$to <= max(coo), ]

    }


    df$sim <- sapply(1:nrow(df), function (i) {

        idx <- df$from[i]:df$to[i]

        sum(apply(ali[, idx], 2, function (v) v[1] == v[2])) / df$len[i] * 100

    })


    return(df)

})

names(L) <- other_seq




# -------------------- PALETTE --------------------
vOTU_repres <- unique(t2$cl_repres[t2$genus == 'delta27'])
vOTU_repres <- vOTU_repres[!(vOTU_repres %in% c('NL_crAss001281', 'NL_crAss001369'))]   # exclude very short contigs

PALETTE <- setNames(
    c('#FFCCCC', '#CCEEFF', '#EEDD88', '#77AADD', '#CC3311', '#EE7733', '#999933', '#AA4499', '#004488', '#009988'),
    c('OLNZ01000083', 'OLNE01000081', 'MT774406.1', 'OHAY01000065', 'NL_crAss001078', 'NL_crAss000239', 'cs_ms_40_1', 'NL_crAss001084', 'NL_crAss000719', 'OBCI01000009')
)
PALETTE <- PALETTE[vOTU_repres]

saveRDS(PALETTE, file = 'delta27_vOTUs_palette.rds')




# -------------------- PLOT --------------------
pdf('RNAp_recombination_plot.pdf', height = 8.75, width = 7.5)

layout(matrix(1:5, nrow = 5), height = c(2,1.8,2,2,1))




# -------------------- plot read coverage --------------------
MAR <- par()$mar
MGP <- par()$mgp
par(mar = MAR + c(-3,2,-1,-1), las = 1, mgp = MGP + 1)

plot(
    NA,
    xlim = c(0, 110000),
    ylim = c(1, 10^3),
    xlab = '',
    ylab = 'Mean read depth + 1',
    cex.lab = 1.2,
    xaxt = 'n', xaxs = 'i',
    yaxt = 'n', log = 'y',
    bty = 'n'
)


par(mgp = MGP)
axis(side = 1, line = 0.5, at = seq(0, 110000, 10000), labels = seq(0, 110, 10), xpd = TRUE)
par(mgp = MGP + 1)
axis(side = 2, at = 10^c(0:3), labels = as.expression(sapply(0:3, function (i) bquote(10^.(i)))))


invisible(
    lapply(R[['OLNE01000081']], function (v) {

        lines(
            x = w_center,
            y = v + 1,
            col = adjustcolor('grey25', alpha.f = 0.2)
        )
}))


mtext('A', side = 3, line = 1, adj = -0.12, family = 'sans', cex = 2)




# -------------------- plot genome map --------------------
par(mar = MAR + c(-2.5,2,0.8,-1))

plot_genome_map(tab, dom, COL, len = t1$cl_member_len[t1$cl_member == 'OLNE01000081'], Xmax = 110000)


text(
    labels = c('capsid', 'transcription', 'replication'),
    x = c(17, 58.8, 80) * 1000,
    y = 7.5,
    font = 2,
    cex = 1.2,
    xpd = TRUE
)


D1 <- c('TerL', 'portal', 'MCP', 'Tstab', 'RNAP', 'VP02740', 'primase', 'PolB')
D2 <- c('TerL', 'portal', 'MCP', 'Tstab', 'RNApBprime', 'VP02740', 'primase', 'DNApB')
text(
    labels = D1,
    col = COL[D2],
    x = c(7, 14.5, 22, 43.8, 58.8, 69.9, 81, 90) * 1000,
    y = -0.5,
    font = 2,
    cex = 1.2,
    xpd = TRUE
)


mtext('B', side = 3, line = 2, adj = -0.12, family = 'sans', cex = 2)




# -------------------- plot nucleotide identity 1 --------------------
# OLNE01000081 vOTU
par(mar = MAR + c(-3,2,-1,-1))

plot(
    NA,
    xlim = c(0, 110000),
    ylim = c(0, 100),
    xlab = '',
    ylab = 'Nucleotide identity, %',
    cex.lab = 1.2,
    xaxt = 'n',
    xaxs = 'i',
    bty = 'n'
)


par(mgp = MGP)
axis(side = 1, line = 0.5, at = seq(0, 110000, 10000), labels = seq(0, 110, 10), xpd = TRUE)


SELE <- t1$cl_member[t1$cl_repres == 'OLNE01000081']
SELE <- SELE[SELE %in% names(L)]


invisible(
    lapply(SELE, function (x) {

        df <- L[[x]]

        lines(
            x = df$at,
            y = df$sim,
            lwd = 0.5,
            col = PALETTE[ t1$cl_repres[t1$cl_member == x] ]
        )
    })
)


mtext('C', side = 3, line = 1, adj = -0.12, family = 'sans', cex = 2)




# -------------------- plot nucleotide identity 2 --------------------
# other delta27 vOTUs
par(mgp = MGP + 1)

plot(
    NA,
    xlim = c(0, 110000),
    ylim = c(0, 100),
    xlab = '',
    ylab = 'Nucleotide identity, %',
    cex.lab = 1.2,
    xaxt = 'n',
    xaxs = 'i',
    bty = 'n'
)


par(mgp = MGP)
axis(side = 1, line = 0.5, at = seq(0, 110000, 10000), labels = seq(0, 110, 10), xpd = TRUE)
mtext('OLNE01000081, kb', side = 1, line = 3.25, xpd = TRUE, cex = 0.8)


SELE <- unlist(lapply(
    c('OLNZ01000083', 'MT774406.1', 'OHAY01000065', 'NL_crAss000239', 'cs_ms_40_1', 'NL_crAss001084', 'NL_crAss000719', 'OBCI01000009', 'NL_crAss001078'),
    function (x) t1$cl_member[t1$cl_repres == x]
))
SELE <- SELE[SELE %in% names(L)]


invisible(
    lapply(SELE, function (x) {

        df <- L[[x]]

        lines(
            x = df$at,
            y = df$sim,
            lwd = 0.5,
            col = PALETTE[ t1$cl_repres[t1$cl_member == x] ]
        )
    })
)


mtext('D', side = 3, line = 1, adj = -0.12, family = 'sans', cex = 2)




# -------------------- plot legend --------------------
par(mar = c(1,1,3,1))

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)

legend(
    'center',
    lty = 1,
    lwd = 1.5,
    col = PALETTE,
    legend = paste(names(PALETTE), 'vOTU'),
    ncol = 5,
    cex = 0.9
)

dev.off()
