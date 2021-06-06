# depends on R packages ape & latex2exp


plot_genera <- function (
    tr,
    DF,
    X1,
    X2,
    sele_labels = c(),
    sele_labels_2shift = c(),
    all_labels = FALSE,
    cex_labels = 1,
    pos_labels = NA) {

    DF <- DF[tr$tip.label, ]


    rle <- rle(DF$genus)


    ### grey bars
    for (i in seq_along(rle$values)) {

        if (rle$values[i] == '') { next }

        from <- ifelse(i==1, 1, sum(rle$lengths[1:(i-1)])+1)
        to <- from + rle$lengths[i] - 1

        rect(
            xleft = X1,
            xright = X2,
            ybottom = from - 0.5,
            ytop = to + 0.5,
            col = ifelse(i%%2 == 0, 'grey80', 'grey90'),
            border = NA
        )

    }


    ### all labels
    if (all_labels) {

        for (i in seq_along(rle$values)) {

            if (rle$values[i] == '') { next }

            from <- ifelse(i==1, 1, sum(rle$lengths[1:(i-1)])+1)
            to <- from + rle$lengths[i] - 1

            g <- sub('[0-9]+$', '', rle$values[i])
            n <- sub('^[a-z]+', '', rle$values[i])

            text(
                x = ifelse(is.na(pos_labels), (X1 + X2) / 2, pos_labels),
                y = (from + to) / 2,
                labels = TeX(paste0('$\\', g, '$', n)),
                cex = cex_labels
            )

        }

    }


    ### selected labels
    for (s in sele_labels) {

        idx <- which(rle$values == s)

        g <- sub('[0-9]+$', '', s)
        n <- sub('^[a-z]+', '', s)

        X0 <- X1 - (X2 - X1) * ifelse(s %in% sele_labels_2shift, 1.1, 0.2)
    
        V <- sapply(idx, function (i) sum(rle$lengths[1:(i-1)]) + rle$lengths[i]/2)
        Y0 <- (min(V) + max(V)) / 2

        text(
            x = X0,
            y = Y0,
            adj = 1.25,
            labels = TeX(paste0('$\\', g, '$', n)),
            cex = cex_labels
        )

        if (length(V) > 1) {

            for (y in V) { lines(x = c(X0, X1), y = c(Y0, y), lwd = 0.5) }

        }

    }

}
