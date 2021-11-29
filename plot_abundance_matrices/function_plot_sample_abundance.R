plot_sample_abundance <- function (i, j, matrix, pal) {

    if (matrix[i, j] != 0) {

        rect(
            xleft = j - 0.5,
            xright = j + 0.5,
            ybottom = i - 0.5,
            ytop = i + 0.5,
            col = pal$col[pal$Above < matrix[i, j] & pal$BelowEqual >= matrix[i, j]],
            border = NA
        )
    }
}


plot_sample_abundance <- Vectorize(plot_sample_abundance, vectorize.args = c('i', 'j'))
