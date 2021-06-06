.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
sessionInfo()



### input parameters
option_list = list(make_option('--x_axis_log', type = 'logical'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
M1 <- read.table(
    '/data/umcg-tifn/crAss_analysis/map_reads/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

M1 <- M1[!grepl('\\|s__', rownames(M1)), ]


M2 <- read.table(
    '/data/umcg-tifn/MetaPhlAn_4cohorts/LLD_LLD2_300OB_IBD_merged_abundance_table.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

colnames(M2) <- sub('_metaphlan$', '', colnames(M2))


t <- read.table(
    '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


z <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/samples_ibd_no_stoma.txt',
    header = FALSE,
    stringsAsFactors = FALSE
)[, 1]



### prepare data
L1 <- list()
L2 <- list()

for (k in unique(t$cohort)) {

    samples <- t$sample_id[ t$cohort == k ]
    if (k == 'IBD') {
        samples <- samples[samples %in% z]
        cat('\n', length(samples), 'IBD samples from subjects w/o stoma were included.\n\n')
    }
    samples <- sapply(samples, function (x) ifelse(grepl('^[0-9]', x), paste0('X', x), x))


    L1[[k]] <- M1[ grep('\\|g__delta27$', rownames(M1)), samples ]

    L2[[k]] <- M2[ grep('\\|s__Prevotella_copri$', rownames(M2)), samples ]

}


COL <- c(LLD = 'forestgreen', LLD2 = 'lightcoral', OB = 'royalblue1', IBD = 'goldenrod1')
COL <- sapply(COL, function (x) adjustcolor(x, alpha.f = 0.5))


if (opt$x_axis_log) {
    Xmin <- 10^(-6)
    Xmax <- 10^ceiling(log10(max(unlist(L1))))
} else {
    Xmin <- 0
    Xmax <- max(pretty(unlist(L1)))
}

Ymax <- max(pretty(unlist(L2)))



### plot data
pdf(
    paste0('delta27_host_corr', ifelse(opt$x_axis_log, '_LOG', ''), '.pdf'),
    width = 3.3,
    height = 3.3
)

par(mar = par()$mar + c(-1.5,-0.5,-3,-1), mgp = par()$mgp + c(-1,-0.5,0), ps = 10)

plot(
    NA,
    xlim = c(Xmin, Xmax),
    ylim = c(0, Ymax),
    xlab = ifelse(
        opt$x_axis_log,
        expression(paste(delta, '27 crAss-like phages abundance + ', 10^-6)),
        expression(paste(delta, '27 crAss-like phages abundance'))
    ),
    ylab = expression(paste(italic('Prevotella copri'), ' abundance, %')),
    tck = -0.03, las = 1,
    xaxt = ifelse(opt$x_axis_log, 'n', 's'),
    log = ifelse(opt$x_axis_log, 'x', '')
)


for (k in unique(t$cohort)) {

    points(
        x = L1[[k]] + Xmin,
        y = L2[[k]],
        pch = 21,
        bg = COL[sub('300OB', 'OB', k)],
        col = COL[sub('300OB', 'OB', k)],
        cex = 0.5
    )

}


if (opt$x_axis_log) {

    N <- seq(from = -6, to = 0, by = 2)
    axis(side = 1, at = 10^N, labels = as.expression(sapply(N, function (i) bquote(10^.(i)))))

}


legend(
    'topright',
    pch = 21,
    cex = 0.5,
    col = COL,
    pt.bg = COL,
    legend = c('LLD', 'LLD follow-up', '300OB', 'IBD')
)

dev.off()
