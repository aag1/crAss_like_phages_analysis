plot_genome_map <- function (tab, dom, col, len, Xmax) {

    ### canvas
    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(6.5, 0.5),
        xaxs = 'i', yaxs = 'i', 
        axes = FALSE, ann = FALSE, bty = 'n'
    )



 	### ORFs
	for (i in 1:nrow(tab)) {

		from   <- tab[i, 2]
		to     <- tab[i, 3]
        strand <- tab[i, 4]


        if (strand == 1) {
            frame <- ((from - 1) %% 3) + 1
            label <- paste0('+', frame)
            Y <- frame
        }
        if (strand == -1) {
            frame <- ((from - 1) %% 3) + 1
            label <- paste0('-', frame)
            Y <- frame + 3
        }


        # grey ORF rectangle
		rect(
            xleft = from,
            xright = to,
            ybottom = Y + 0.5,
            ytop = Y - 0.5,
            col = 'grey98',
            border = NA
        )


        # protein domains
        idx <- which(rownames(dom) == tab[i, 1])
        if (length(idx) == 1) {
            for (d in colnames(dom)) {

                if (dom[idx, d] == '-') { next }

                S <- strsplit(dom[idx, d], ';')[[1]]
                L <- lapply(strsplit(S, '-'), as.numeric)
                d1 <- unlist(lapply(L, function (v) v[1]))
                d2 <- unlist(lapply(L, function (v) v[2]))

                if (strand == 1) {
                    d_from <- (from - 1) + d1*3 - 2
                    d_to <- (from - 1) + d2*3
                }
                if (strand == -1) {
                    d_from <- (to + 1) - d2*3
                    d_to <- (to + 1) - d1*3 + 2
                }

                rect(
                    xleft = d_from,
                    xright = d_to,
                    ybottom = Y + 0.5,
                    ytop = Y - 0.5,
                    col = col[d],
                    border = NA
                )
            }
        }


        # black ORF border
		rect(
            xleft = from,
            xright = to,
            ybottom = Y + 0.5,
            ytop = Y - 0.5,
            col = NA,
            border = 'black'
        )

    }



	### genome
	rect(
        xleft = 0,
        xright = len, 
        ybottom = 6.5,
        ytop = 0.5,
        col = NA,
        border = 'black'
    )

	text(
        x = -0.01 * Xmax,
        y = 1:6,
        labels = c('+1', '+2', '+3', '-1', '-2', '-3'),
        adj = 1,
        xpd = TRUE
    ) 

}
