.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(latex2exp)
sessionInfo()
source('/home/umcg-agulyaeva/crAss_analysis/genome_maps/function_plot_genome_map.R')



### input parameters
option_list = list(
	make_option('--domainsF'),
	make_option('--genomesF'),
    make_option('--GC_skewF'),
    make_option('--CGC_skewF'),
    make_option('--AT_skewF'),
    make_option('--CAT_skewF'),
    make_option('--groupsF'),
	make_option('--outF'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
d1 <- read.table(
    opt$domainsF,
    sep = '\t',
    header = FALSE,
    comment.char = '',
    stringsAsFactors = FALSE)

COL <- setNames(d1[, 2], d1[, 1])


d2 <- read.table(
    opt$genomesF,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE)

G <- unique(d2$cl_repres)


skew <- list(
    GC = readRDS(opt$GC_skewF),
    AT = readRDS(opt$AT_skewF)
)

c_skew <- list(
    GC = readRDS(opt$CGC_skewF),
    AT = readRDS(opt$CAT_skewF)
)


if (!is.null(opt$groupsF)) {

    d3 <- read.table(
        opt$groupsF,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE)

    G <- rownames(d3)[rownames(d3) %in% G]

}



### plot
pdf(opt$outF, height = 7.5, width = 11.7)


DOM_NAMES <- names(COL)
DOM_NAMES[ DOM_NAMES == 'DNApB' ] <- 'PolB'
DOM_NAMES[ DOM_NAMES == 'RNApBprime' ] <- 'RNApB\''
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
legend('center', fill = COL, legend = DOM_NAMES, ncol = 3)


par(oma = c(2, 1, 3, 1))
layout(matrix(1:5, nrow = 5))


for (contig_name in G) {

    par(mar = c(1, 4, 3, 4))
    for (m in c('c11', 'TAG|q', 'TGA|w')) {

        l <- d2$cl_member_len[ d2$cl_member == contig_name ]

        tabF <- paste0('ANNOTATIONS/', contig_name, '/', contig_name, '_', sub('\\|', '', m), '_AA_coordinates.txt')
        tab <- read.table(tabF, sep = '\t', header = FALSE, stringsAsFactors = FALSE)

        domF <- paste0('ANNOTATIONS/', contig_name, '/', 'crAss_key_doms_in_', contig_name, '_', sub('\\|', '', m), '.txt')
        dom <- read.table(domF, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)


        # genome map
        plot_genome_map(tab, dom, COL, len = l, Xmax = 200000)


        # translation
        FONT <- 1
        if (!is.null(opt$groupsF)) {
            FONT <- ifelse(d3[contig_name, 'translation'] == m, 2, 1)
        }

        text(
            x = 0,
            y = -0.5,
            labels = paste(contig_name, ifelse(m == 'c11', '', m)),
            font = FONT,
            adj = 0,
            cex = 1.5,
            xpd = TRUE
        )


        # axis
        if (m == 'c11') {
            axis(side = 3, outer = TRUE, col = 'grey75', col.axis = 'grey75')
        }


        # group label
        if (m == 'c11' & !is.null(opt$groupsF)) {

            genus <- d3[contig_name, 'genus']

            if (genus != '') {

                gr <- sub('[0-9]+$', '', genus)
                nm <- sub('^[a-z]+', '', genus)

                text(
                    x = par()$usr[2],
                    y = par()$usr[4],
                    labels = TeX(paste0('$\\', gr, '$', nm)),
                    xpd = TRUE,
                    cex = 2.5
                )

            }

        }

    }


    # GC & AT skew
    par(mar = c(2, 4, 2, 4))

    w_center <- seq(from = 501, to = l - 500, by = 200)

    for (x in c('GC', 'AT')) {

        r <- range(c(skew[['GC']][[ contig_name ]], skew[['AT']][[ contig_name ]]), na.rm = TRUE)
        plot(
            x = w_center,
            y = skew[[x]][[ contig_name ]],
            type = 'l',
            col = 'blue',
            xlim = c(0, 200000),
            ylim = c(floor(r[1]/0.2), ceiling(r[2]/0.2))*0.2,
            xaxs = 'i', yaxs = 'i',
            axes = FALSE, ann = FALSE, bty = 'n'
        )
        axis(side = 2, line = par()$mgp[3], col = 'blue', col.axis = 'blue', las = 1)
        mtext(paste0(x, ' skew'), line = par()$mgp[1], side = 2, col = 'blue')

        par(new = TRUE)

        r <- range(c(c_skew[['GC']][[ contig_name ]], c_skew[['AT']][[ contig_name ]]), na.rm = TRUE)
        plot(
            x = w_center,
            y = c_skew[[x]][[ contig_name ]],
            type = 'l',
            col = 'red',
            xlim = c(0, 200000),
            ylim = c(floor(r[1]/20), ceiling(r[2]/20))*20,
            xaxs = 'i', yaxs = 'i',
            axes = FALSE, ann = FALSE, bty = 'n'
        )
        axis(side = 1, line = 1, col = 'grey75', col.axis = 'grey75', labels = (x == 'AT'))
        axis(side = 4, line = par()$mgp[3], col = 'red', col.axis = 'red', las = 1)
        mtext(paste0('C', x, ' skew'), line = par()$mgp[1], side = 4, col = 'red')

    }

}

dev.off()
