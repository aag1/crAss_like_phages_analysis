### are circular crAss contigs enriched in any family?

d <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/NL_crAss_contigs_2.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


k <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


t <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


d$cl_repres <- sapply(d$contig_name, function (x) k$cl_repres[k$cl_member == x])


d$group <- sapply(d$cl_repres, function (x) t$group[t$cl_repres == x])


for (x in c('alpha_gamma', 'beta', 'delta', 'epsilon', 'zeta')) {

    S <- mean(d$length[d$group == x])

    A <- sum(d$group == x)

    B <- sum((d$group == x) & (d$ends_overlap == 1))

    N <- round(B / A * 100, 1)

    cat(x, ':', A, 'contigs,', N, '% circular contigs,', S, 'nt mean contig length\n')

}
#alpha_gamma : 682 contigs, 9.7 % circular contigs, 41944.29 nt mean contig length
#beta : 68 contigs, 27.9 % circular contigs, 77122.57 nt mean contig length
#delta : 400 contigs, 34.5 % circular contigs, 70474.06 nt mean contig length
#epsilon : 335 contigs, 1.5 % circular contigs, 57833.82 nt mean contig length
#zeta : 71 contigs, 38 % circular contigs, 141892.6 nt mean contig length




### correlation between the number of clean reads and the number of (1) crAss contigs, (2) crAss vOTUs per sample
L <- readRDS('/data/umcg-tifn/crAss_analysis/detection_summary/reads_contigs_number.rds')


DF <- NULL
for (x in names(L)) {

    df <- L[[x]]

    DF <- rbind(DF, df)

}


OBJ <- cor.test(x = DF$clean_reads, y = DF$num_crAss_contigs, method = 'spearman')
OBJ
#S = 1638892642, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.1822387


OBJ$p.value
#1.473133e-18



OBJ <- cor.test(x = DF$clean_reads, y = DF$num_seq_cov10pct, method = 'spearman')
OBJ
#S = 1577121307, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.2130608

OBJ$p.value
#6.279732e-25




### most abundant bacterial phyla in 3 cohorts
M <- read.table(
    '/data/umcg-tifn/MetaPhlAn_4cohorts/LLD_LLD2_300OB_IBD_merged_abundance_table.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

colnames(M) <- sub('_metaphlan$', '', colnames(M))


t <- read.table(
    '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


z <- read.table(
    '/home/umcg-agulyaeva/crAss_analysis/crAss_hosts/selected_ibd_samples.txt',
    header = FALSE,
    stringsAsFactors = FALSE
)[, 1]


sele <- grep('^k__Bacteria\\|p__[^\\|]+$', rownames(M), value = TRUE)


for (k in c('LLD', '300OB', 'IBD')) {

    samples <- t$sample_id[ t$cohort == k ]
    if (k == 'IBD') { samples <- samples[samples %in% z] }
    samples <- sapply(samples, function (x) ifelse(grepl('^[0-9]', x), paste0('X', x), x))


    m <- M[sele, samples]


    v <- apply(m, 1, mean)
    v <- sort(v, decreasing = TRUE)
    v <- round(v)


    df <- data.frame(phylum = sub('^k__Bacteria\\|p__', '', names(v)), abundance = v)
    rownames(df) <- NULL


    cat('\n\n', k, '\n')
    print(df[1:5,])

}


# LLD 
#           phylum abundance
#1      Firmicutes        64
#2  Actinobacteria        21
#3   Bacteroidetes        11
#4 Verrucomicrobia         2
#5  Proteobacteria         1


# 300OB 
#           phylum abundance
#1      Firmicutes        61
#2   Bacteroidetes        24
#3  Actinobacteria        11
#4  Proteobacteria         2
#5 Verrucomicrobia         1


# IBD 
#           phylum abundance
#1      Firmicutes        61
#2  Actinobacteria        18
#3   Bacteroidetes        18
#4  Proteobacteria         1
#5 Verrucomicrobia         1
