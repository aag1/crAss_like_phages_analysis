.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(latex2exp)
sessionInfo()



### input parameters
option_list = list(
    make_option('--LLD_depth'),
	make_option('--LLD2_depth'),
	make_option('--OB_depth'),
    make_option('--IBD_depth'),
    make_option('--nt_content'),
    make_option('--groups_tab'),
    make_option('--outF'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read data
L <- list(
    'LLD'   = readRDS(opt$LLD_depth),
    'LLD2'  = readRDS(opt$LLD2_depth),
    '300OB' = readRDS(opt$OB_depth),
    'IBD'  = readRDS(opt$IBD_depth)
)


D <- readRDS(opt$nt_content)


G <- names(D)

sele <- unique(unlist(lapply(L, function (l) names(l))))
G <- G[G %in% sele]

if (!is.null(opt$groups_tab)) {

    tab <- read.table(
        opt$groups_tab,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE)

    G <- rownames(tab)[rownames(tab) %in% G]

}



### plot data
pdf(opt$outF, height = 11.7, width = 8.3)

layout(matrix(1:5, nrow = 5))

par(las = 1, ps = 14)

MAR <- par()$mar

for (g in G) {

    ### window 1001 nt, step 200 nt
    Xcoo <- 501 + 0:(length(D[[g]][[1]]) - 1) * 200


    ### plot depth
    for (x in names(L)) {

        if (x == 'LLD') {
            par(mar = MAR + c(-2.5, 3, -3, 0))
        } else {
            par(mar = MAR + c(-2.5, 3, -3, 0))
        }


        plot(
            NA,
            xlim = c(0, 200000),
            ylim = c(0, 4),
            xlab = paste0(g, ', nt'),
            ylab = paste(ifelse(x == 'LLD2', 'LLD follow-up', x), '\nlog10(1 + mean depth)'))


        if (g %in% names(L[[x]])) {

            invisible(lapply(L[[x]][[g]], function (v) {

                lines(
                    x = Xcoo,
                    y = log10( v + 1 ),
                    col = adjustcolor('grey25', alpha.f = 0.2))
            }))
        }


        if (x == 'LLD' & !is.null(opt$groups_tab)) {

            genus <- tab[g, 'genus']

            if (genus != '') {

                gr <- sub('[0-9]+$', '', genus)
                nm <- sub('^[a-z]+', '', genus)

                text(
                    x = par()$usr[2],
                    y = par()$usr[4],
                    labels = TeX(paste0('$\\', gr, '$', nm)),
                    adj = c(2, 2),
                    cex = 2)
            }
        }
    }


    ### plot nt content
    par(mar = MAR + c(0, 3, -3, 0))

    plot(
        NA,
        xlim = c(0, 200000),
        ylim = c(0, 100),
        xlab = paste0(g, ', nt'),
        ylab = 'nt content, %'
    )

    lines(x = Xcoo, y = D[[g]][['A']], col = 'goldenrod')
    lines(x = Xcoo, y = D[[g]][['T']], col = 'green')
    lines(x = Xcoo, y = D[[g]][['G']], col = 'blue')
    lines(x = Xcoo, y = D[[g]][['C']], col = 'red')

    legend(
        'topright',
        lty = 1,
        col = c('goldenrod', 'green', 'blue', 'red'),
        legend = c('A', 'T', 'G', 'C')
    )

}

dev.off()
