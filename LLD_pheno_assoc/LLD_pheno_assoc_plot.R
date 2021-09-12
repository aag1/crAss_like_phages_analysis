.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(latex2exp)
sessionInfo()




# -------------------- read data --------------------
crass = read.table('../from_Peregrine/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt', header = T, as.is = T)



linkage_file = read.table('../../LLD_LLD2_Info/LLD_GTMicrob.txt', as.is = T)

linkage_file$V2 = paste0('X', linkage_file$V2)

shared_samples = intersect(linkage_file$V2, colnames(crass))

linkage_file2 = linkage_file[match(shared_samples, linkage_file$V2), ]
# 'match' returns a vector of the positions of (first) matches of its first argument in its second



crass1 = crass[grep('\\|s__', rownames(crass), invert = T), ]

crass2 = t(crass1[, match(shared_samples, colnames(crass1))])

crass3 = crass2[, colSums(crass2 > 0) > nrow(crass2) * 0.05]



pheno <- data.frame(NULL)

pheno_files <- list.files(path = '../../LLD_LLD2_Info/Pheno_science_imputed_1135/', full.names = T)

for (f in pheno_files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != 'LLDEEPid']

    if (ncol(pheno) == 0) { pheno <- df } else { pheno <- cbind(pheno, df) }

}

pheno2 = pheno[match(linkage_file2$V1, rownames(pheno)), ]



tab <- read.table(
    'LLD_pheno_assoc.txt',
    sep = '\t',
    header = T,
    stringsAsFactors = F
)

tab <- tab[tab$FDR.bin < 0.05, ]




# -------------------- functions --------------------
stats_lab <- function (PH, TX, tab) {

    TX1 <- sub('^o__', '', TX)
    TX1 <- gsub('\\|[a-z]__', '|', TX1)
    idx <- which(tab$pheno == PH & tab$taxon == TX1)

    r <- tab[idx, 'Rsp.bin']
    r <- format(round(r, 3), nsmall = 3)

    p <- tab[idx, 'P.bin']
    p <- formatC(p, format = 'e', digits = 2)

    FDR <- tab[idx, 'FDR.bin']
    FDR <- formatC(FDR, format = 'e', digits = 2)

    lab <- paste0('r = ', r, '\nP-value = ', p, '\nFDR = ', FDR)

    return(lab)

}



taxo_lab <- function (TX) {

    TX2 <- strsplit(TX, '\\|')[[1]]
    TX2 <- TX2[ length(TX2) ]
    TX2 <- sub('^[a-z]__', '', TX2)

    if (TX2 == 'alpha_gamma') {

        lab <- '$\\alpha\\gamma$'

    } else {

        g <- sub('[0-9]+$', '', TX2)
        n <- sub('^[a-z_]+', '', TX2)
        lab <- paste0('$\\', g, '$', n)

    }

    return(lab)

}




# -------------------- plot data --------------------
pdf('LLD_pheno_assoc_plot.pdf', height = 7.5, width = 7.5)

layout(matrix(c(rep(1, 40), rep(2, 10), rep(3, 20), rep(4, 10), rep(5, 7), rep(6, 5), rep(7, 4), rep(8, 8), rep(9, 8), rep(10, 8)), nrow = 3, ncol = 40, byrow = TRUE))

MAR0 <- par()$mar



### continuous phenotypes
for (PH in c('BioMK_ChromograninA_log', 'BioMK_BetaDefensin2_log', 'meat_log', 'carbohydrates.total_log')) {

    if (PH == 'BioMK_ChromograninA_log') {

        taxo <- c('o__crAss', 'o__crAss|f__alpha_gamma', 'o__crAss|f__beta', 'o__crAss|f__alpha_gamma|g__gamma1', 'o__crAss|f__alpha_gamma|g__gamma4')
        YLIM <- c(0, 3)
        YLAB <- 'CgA, nmol/g'
        MAR <- MAR0 + c(0, 1, 0, -1)

    }

    if (PH == 'BioMK_BetaDefensin2_log') {

        taxo <- 'o__crAss|f__alpha_gamma|g__gamma1'
        YLIM <- c(0, 3)
        YLAB <- 'HBD-2, ng/g'
        MAR <- MAR0 + c(0, 1, 0, -1)

    }

    if (PH == 'meat_log') {

        taxo <- c('o__crAss|f__epsilon', 'o__crAss|f__epsilon|g__epsilon1')
        YLIM <- c(-1, 3)
        YLAB <- 'Meat, g/day'
        MAR <- MAR0 + c(0, 3, 0, 0)

    }

    if (PH == 'carbohydrates.total_log') {

        taxo <- 'o__crAss|f__alpha_gamma|g__alpha20'
        YLIM <- c(1, 3)
        YLAB <- 'Carbohydrates, g/day'
        MAR <- MAR0 + c(0, 1, 0, -1)

    }


    par(mar = MAR)

    plot(
        NA,
        xlim = c(0.5, 0.5 + length(taxo) * 1.25),
        ylim = YLIM,
        xaxs = 'i',
        yaxs = 'i',
        xaxt = 'n',
        yaxt = 'n',
        bty = 'n',
        ann = FALSE
    )

    I <- seq(from = YLIM[1], to = YLIM[2], by = 1)
    axis(
        side = 2,
        at = I,
        labels = as.expression(sapply(I, function (i) bquote(10^.(i)))),
        las = 1
    )

    mtext(
        text = YLAB,
        side = 2,
        line = 3,
        cex = 0.8
    )


    X <- 1

    for (TX in taxo) {

        V1 <- pheno2[crass3[, TX] == 0, PH]
        V2 <- pheno2[crass3[, TX] > 0, PH]


        boxplot(
            V1,
            pch = 20,
            add = TRUE,
            at = X,
            yaxt = 'n',
            frame = FALSE
        )

        stripchart(
            V1,
            method = 'jitter',
            jitter = 0.2,
            pch = 20,
            cex = 0.01,
            col = adjustcolor('forestgreen', alpha.f = 0.5),
            vertical = TRUE,
            add = TRUE,
            at = X
        )


        boxplot(
            V2,
            pch = 20,
            add = TRUE,
            at = X + 0.5,
            yaxt = 'n',
            frame = FALSE
        )

        stripchart(
            V2,
            method = 'jitter',
            jitter = 0.2,
            pch = 20,
            cex = 0.01,
            col = adjustcolor('forestgreen', alpha.f = 0.5),
            vertical = TRUE,
            add = TRUE,
            at = X + 0.5
        )


        mtext(
            text = stats_lab(PH, TX, tab),
            side = 3,
            line = 1,
            at = X + 0.25,
            cex = 0.6
        )

        mtext(
            text = c('-', '+'),
            font = 2,
            side = 1,
            line = 0.5,
            at = X + c(0, 0.5)
        )

        mtext(
            text = ifelse(TX == 'o__crAss', 'all', TeX(taxo_lab(TX))),
            side = 1,
            line = 2,
            at = X + 0.25
        )


        X <- X + 1.25

    }

}



### categorical phenotypes
PH <- 'how_often_coffee'

colfunc <- colorRampPalette(c('grey25', 'grey95'))

for (TX in c('o__crAss', 'o__crAss|f__beta')) {

    if (TX == 'o__crAss') {

        par(mar = MAR0 + c(0, 1, 0, -2), mgp = c(3, 1, 0.1))

    } else {

        par(mar = MAR0 + c(0, -1.8, 0, -2))

    }


    M <- matrix(0, nrow = 7, ncol = 2)

    for (i in 0:6) {

        M[i + 1, 1] <- sum(pheno2[crass3[, TX] == 0, PH] == i)

        M[i + 1, 2] <- sum(pheno2[crass3[, TX] > 0, PH] == i)

    }

    M <- apply(M, 2, function (v) v / sum(v) * 100)

    bp <- barplot(
        M,
        col = colfunc(7),
        axes = (TX == 'o__crAss'),
        las = 1
    )


    if (TX == 'o__crAss') {

       mtext(
           text = '% samples',
           side = 2,
           line = 3,
           cex = 0.8
       )

    }

    mtext(
        text = stats_lab(PH, TX, tab),
        side = 3,
        line = 1,
        at = sum(bp) / 2,
        cex = 0.6
    )

    mtext(
        text = c('-', '+'),
        font = 2,
        side = 1,
        line = 0.5,
        at = bp
    )

    mtext(
        text = ifelse(TX == 'o__crAss', 'all', TeX(taxo_lab(TX))),
        side = 1,
        line = 2,
        at = sum(bp) / 2
    )

}


par(mar = c(MAR0 + c(0, -3, 0, -1)))

plot(
    NULL,
    xlim = c(0, 1),
    ylim = c(0, 1),
    axes = FALSE,
    ann = FALSE
)

legend(
    'center',
    fill = rev(colfunc(7)),
    legend = 6:0,
    title = 'coffee\nintake\nfrequency\n',
    bty = 'n',
    xpd = TRUE
)



### binary phenotypes
par(mar = MAR0 + c(0, 2, 0, 0))

TX <- 'o__crAss|f__alpha_gamma'

for (PH in c('IBS', 'stomach_ulcer', 'laxatives')) {

    if (PH == 'IBS') {

        YLIM <- c(0, 20)
        YLAB <- 'IBS, % samples'

    }

    if (PH == 'stomach_ulcer') {

        YLIM <- c(0, 10)
        YLAB <- 'Stomach ulcer, % samples'

    }

    if (PH == 'laxatives') {

        YLIM <- c(0, 10)
        YLAB <- 'Laxatives usage, % samples'

    }


    n1 <- sum(pheno2[crass3[, TX] == 0, PH] == 1) / sum(crass3[, TX] == 0) * 100
    n2 <- sum(pheno2[crass3[, TX] > 0, PH] == 1) / sum(crass3[, TX] > 0) * 100

    bp <- barplot(
        c(n1, n2),
        ylim = YLIM,
        yaxs = 'i',
        las = 1
    )


    mtext(
        text = YLAB,
        side = 2,
        line = 2.5,
        cex = 0.8
    )

    mtext(
        text = stats_lab(PH, TX, tab),
        side = 3,
        line = 1,
        at = sum(bp) / 2,
        cex = 0.6
    )

    mtext(
        text = c('-', '+'),
        font = 2,
        side = 1,
        line = 0.5,
        at = bp
    )

    mtext(
        text = ifelse(TX == 'o__crAss', 'all', TeX(taxo_lab(TX))),
        side = 1,
        line = 2,
        at = sum(bp) / 2
    )

}

dev.off()
