# depends on R packages ape & latex2exp, function_plot_tree_support_pch.R


plot_crass_tree <- function (
    tr,
    DF,
    xMax,
    bootstrap = FALSE,
    tip_labels = FALSE,
    cex_tip_labels = 1,
    cex_group_labels = 1,
    adj_group_labels = list()) {

    DF <- DF[tr$tip.label, ]



    ### color branches by group
    PALETTE <- setNames(
        c('#AA4488', '#4477AA', '#44AA77', '#AAAA44', '#AA7744'),
        c('alpha_gamma', 'beta', 'delta', 'epsilon', 'zeta')
    )

    EDGE_COL <- rep('black', Nedge(tr))

    G <- unique(DF$group)
    G <- G[G != '']

    for (g in G) {

        clade_tips <- rownames(DF)[ DF$group == g ]

        clade_edges <- which.edge(tr, group = clade_tips)

        EDGE_COL[ clade_edges ] <- PALETTE[g]

    }



    ### plot tree
    plot.phylo(
        tr,
        edge.width = 1.5,
        edge.color = EDGE_COL,
        show.tip.label = tip_labels,
        label.offset = 0.01,
        cex = cex_tip_labels,
        x.lim = c(0, xMax)
    )

    add.scale.bar()



    ### bootstrap support
    if (bootstrap) {

        d <- as.data.frame(tr$edge, stringsAsFactors = FALSE)
        colnames(d) <- c('nodeA', 'nodeD')
        d$edge_len <- tr$edge.length
        d <- d[ d$nodeD > Ntip(tr), ]
        d <- d[ d$edge_len > 0.05, ]

        tree_support_pch(
            tree = tr,
            node = d$nodeD,
            bs_pos = 'topleft',
            bs_name = 'BP',
            max = 100
        )

    }



    ### group labels
    for (g in G) {

        clade_tips <- rownames(DF)[ DF$group == g ]

        clade_mrca <- getMRCA(tr, tip = clade_tips)

        latex_label <- ifelse(g == 'alpha_gamma', '$\\alpha\\gamma$', paste0('$\\', g, '$'))

        if (g %in% names(adj_group_labels)) { ADJ <- adj_group_labels[[g]] } else { ADJ <- c(1.5, -1) }

        nodelabels(
            node = clade_mrca,
            col = PALETTE[g],
            text = TeX(latex_label),
            frame = 'none',
            adj = ADJ,
            cex = cex_group_labels
        )

    }

}
