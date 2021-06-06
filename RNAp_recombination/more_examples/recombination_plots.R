.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(bio3d)
library(latex2exp)
source('/home/umcg-agulyaeva/crAss_analysis/genome_maps/function_plot_genome_map.R')
sessionInfo()




# -------------------- read data --------------------
R_lld <- readRDS('/data/umcg-tifn/crAss_analysis/map_reads/LLD_crAss_depth.rds')
R_300OB <- readRDS('/data/umcg-tifn/crAss_analysis/map_reads/300OB_crAss_depth.rds')


d <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/genome_maps/selected_domains.txt',
    sep = '\t',
    header = FALSE,
    comment.char = '',
    stringsAsFactors = FALSE)

d <- d[!(d[, 1] %in% c('BACON', 'RT_G2')), ]

COL <- setNames(d[, 2], d[, 1])


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




# -------------------- PLOT --------------------
pdf('RNAp_recombination_alpha_beta_zeta.pdf', height = 7.5, width = 11.7)

layout(matrix(1:4, nrow = 4), height = c(2,1.8,2,1))

MAR <- par()$mar
MGP <- par()$mgp


GENUS <- c('alpha6', 'beta8', 'zeta9')
GENOME <- c('OLOZ01000098', 'IAS_virus_KJ003983', 'NL_crAss000848')
X_max_kb <- c(110, 110, 180)


for (i in 1:3) {

genus <- GENUS[i]
genome <- GENOME[i]
tr <- sub('\\|', '', t2$translation[t2$cl_repres == genome])


tabF <- paste0('/data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/', genome, '/', genome, '_', tr, '_AA_coordinates.txt')
tab <- read.table(tabF, sep = '\t', header = FALSE, stringsAsFactors = FALSE)


domF <- paste0('/data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/', genome, '/crAss_key_doms_in_', genome, '_', tr, '.txt')
dom <- read.table(domF, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)


MSA <- read.fasta(paste0(genus, '_genomes_circular_MSA.fasta'))$ali




# -------------------- calculations for simplot --------------------
sele_seq <- genome
other_seq <- rownames(MSA)[rownames(MSA) != sele_seq]
if (sele_seq == 'IAS_virus_KJ003983') { other_seq <- 'NL_crAss001291' }


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


    df$sim <- sapply(1:nrow(df), function (i) {

        idx <- df$from[i]:df$to[i]

        sum(apply(ali[, idx], 2, function (v) v[1] == v[2])) / df$len[i] * 100

    })


    return(df)

})

names(L) <- other_seq




# -------------------- PALETTE --------------------
vOTU_repres <- unique(t1$cl_repres[t1$cl_member %in% other_seq])

if (i == 1) { col <- c('indianred1', 'skyblue2', '#004488') }
if (i == 2) { col <- c('#77AADD') }
if (i == 3) { col <- c('#999933', '#AA4499') }

PALETTE <- setNames(col, vOTU_repres)




# -------------------- plot read coverage --------------------
par(mar = MAR + c(-3,2,-2,-1), las = 1, mgp = MGP + 1)

plot(
    NA,
    xlim = c(0, 180000),
    ylim = c(1, 10^3),
    xlab = '',
    ylab = 'Mean read depth + 1',
    cex.lab = 1.2,
    xaxt = 'n', xaxs = 'i',
    yaxt = 'n', log = 'y',
    bty = 'n'
)


par(mgp = MGP)
axis(side = 1, line = 0.5, at = seq(0, X_max_kb[i] * 1000, 20000), labels = seq(0, X_max_kb[i], 20), xpd = TRUE)
par(mgp = MGP + 1)
axis(side = 2, at = 10^c(0:3), labels = as.expression(sapply(0:3, function (i) bquote(10^.(i)))))


R <- R_lld
if (i == 2) { R <- R_300OB }

invisible(
    lapply(R[[genome]], function (v) {

        lines(
            x = w_center,
            y = v + 1,
            col = adjustcolor('grey25', alpha.f = 0.2)
        )
}))


# genus label
gr <- sub('[0-9]+$', '', genus)
nm <- sub('^[a-z]+', '', genus)

text(
    x = par()$usr[2],
    y = 10^par()$usr[4],
    labels = TeX(paste0('$\\', gr, '$', nm)),
    adj = c(1, 0.5),
    xpd = TRUE,
    cex = 2.5
)




# -------------------- plot genome map --------------------
par(mar = MAR + c(-2.5,2,0.8,-1))

plot_genome_map(tab, dom, COL, len = t1$cl_member_len[t1$cl_member == genome], Xmax =180000)

if (i != 3) {

    legend(
        x = 140000,
        y = par()$usr[3] + diff(par()$usr[3:4]) / 2,
        fill = COL,
        legend = names(COL),
        xjust = 0.5,
        yjust = 0.5,
        ncol = 2
    )

}




# -------------------- plot nucleotide identity --------------------
par(mar = MAR + c(-3,2,-2,-1))

plot(
    NA,
    xlim = c(0, 180000),
    ylim = c(0, 100),
    xlab = '',
    ylab = 'Nucleotide identity, %',
    cex.lab = 1.2,
    xaxt = 'n',
    xaxs = 'i',
    bty = 'n'
)


par(mgp = MGP)
axis(side = 1, line = 0.5, at = seq(0, X_max_kb[i] * 1000, 20000), labels = seq(0, X_max_kb[i], 20), xpd = TRUE)
mtext(paste0(genome, ', kb'), side = 1, line = 3.5, at = X_max_kb[i] / 2 * 1000, xpd = TRUE, cex = 0.8)


invisible(
    lapply(names(L), function (x) {

        df <- L[[x]]

        lines(
            x = df$at,
            y = df$sim,
            lwd = 0.5,
            col = PALETTE[ t1$cl_repres[t1$cl_member == x] ]
        )
    })
)




# -------------------- plot legend --------------------
par(mar = c(1,1,3,1))

plot(NULL, xlim = c(0, 180), ylim = c(0, 1), axes = FALSE, ann = FALSE)

legend(
    x = X_max_kb[i] / 2,
    y = 0.5,
    lty = 1,
    lwd = 1.5,
    col = PALETTE,
    legend = paste(names(PALETTE), 'vOTU'),
    xjust = 0.5,
    horiz = TRUE
)

}

dev.off()
