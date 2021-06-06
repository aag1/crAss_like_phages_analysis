# function to plot bootstrap support


tree_support_pch <- function (tree, node, bs_pos, bs_name, max) {

    dat <- data.frame(node, stringsAsFactors = FALSE)
    dat$bs <- sapply(dat$node, function (n) tree$node.label[n - Ntip(tree)])
    dat <- dat[dat$bs != '', ]
    dat$bs <- as.numeric(dat$bs)

	if (max == 1) { dat$bs <- dat$bs * 100 }

	dat$bg <- sapply(dat$bs, function (x) {
					if (x >= 90)          { return('black')  }
					if (x < 90 & x >= 70) { return('grey75') }
					if (x < 70)           { return('white')  }
	})

	nodelabels(node = dat$node, pch = 21, bg = dat$bg, col = 'black')

	if (!is.na(bs_pos)) {
		legend(
            bs_pos,
			legend = eval(substitute(expression('90%' <= bs, '70%' <= bs * ' < 90%', bs < '70%'), list(bs=bs_name))),
			pch = 21,
            pt.bg = c('black', 'grey', 'white'),
			cex = 1,
            xpd = TRUE
        )
	}

}
