sessionInfo()




# ------------------------------ CRISPR ------------------------------ #
tab <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_hosts/crAss_host_pairs.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


tab$host_genus <- unlist(lapply(strsplit(tab$host_taxonomy, ';'), function (v) {

    x <- sub('^uncultured ', '', rev(v)[1])
    y <- strsplit(x, ' ')[[1]][1]
    return(y)

}))


length(unique(tab$host_genus))  # number of host genera
# 10

length(unique(tab$phage_vOTU))  # number of crAss vOTUs
# 96

length(unique(tab$phage_genus)) # number of crAss genus-level clusters
# 32

sum(tab$host_genus %in% c('Bacteroides', 'Prevotella', 'Porphyromonas', 'Parabacteroides')) / nrow(tab) * 100
# 97.79893


idx <- which(tab$phage_genus == 'delta27')
d27_phage_id <- unique(tab$phage_id[idx])
d27_host_id <- unique(tab$host_id[idx])


d27_host_id
# CP042464.1


t1 <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_hosts/crAss_host_pairs.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


d27_spacer <- t1$spacer[(t1$phage %in% d27_phage_id) & (t1$host %in% d27_host_id)]
d27_spacer <- unique(d27_spacer)
d27_spacer
# [1] "spacer286428;spacer089709"             
# [2] "spacer286428"                          
# [3] "spacer111241;spacer286428"             
# [4] "spacer111241"                          
# [5] "spacer275021;spacer089709;spacer286428"


# vOTUs linked to hosts by CRISPR, that are not on the tree
library(ape)

tr <- read.tree('/data/umcg-tifn/crAss_analysis/build_tree/MSA_crAss_TerL_MidpointRooted.2.newick')

tab[!(tab$phage_vOTU %in% tr$tip.label), c(1:3, 7)]
# phage_id phage_vOTU phage_genus host_genus
# phage_6_3  phage_6_3     delta37 Prevotella
# phage_6_3  phage_6_3     delta37 Prevotella
# phage_6_3  phage_6_3     delta37 Prevotella
# phage_6_3  phage_6_3     delta37 Prevotella
# phage_13_2 phage_13_2     delta34 Prevotella
# phage_14_3 phage_14_3      zeta29 Prevotella
# phage_14_3 phage_14_3      zeta29 Prevotella
# phage_14_3 phage_14_3      zeta29 Prevotella
# phage_14_3 phage_14_3      zeta29 Prevotella
# phage_6_2  phage_6_2     delta35 Prevotella




# ------------------------------ Guerin et al. 2018 genera ------------------------------ #
b <- read.csv(
    '/data/umcg-tifn/DATABASES/data_Guerin_2018/Data_S1/crAss_with_yutin_seqs_nodes.csv',
    header = TRUE,
    stringsAsFactors = FALSE
)

k <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

g <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


b <- b[!is.na(b$Candidate_Genus), ]

b$genus_level_cluster <- NA

for (i in 1:nrow(b)) {

    cl_member <- tolower(b$CONTIG_NAME[i])

    idx <- grep(paste0('^',cl_member,'(_1){0,1}$'), k$cl_member, ignore.case = TRUE)
    if (length(idx) == 0) { next }

    cl_repres <- k$cl_repres[idx]

    b$genus_level_cluster[i] <- g$genus[g$cl_repres == cl_repres]

}

nrow(b)
# 249


for (g1 in 1:10) {

    g2 <- unique(b$genus_level_cluster[b$Candidate_Genus == g1])
    g2 <- sort(g2)

    cat(g1, ':', paste(g2, collapse = ', '), '\n\n')

}
# 1 : alpha20 
#
# 2 : gamma1 
#
# 3 : alpha10, alpha13, alpha15, alpha17, alpha18, alpha21, alpha23, alpha24 
#
# 4 : alpha26, alpha27 
#
# 5 : gamma3, gamma4 
#
# 6 : beta10, beta13, beta16, beta2, beta23, beta26, beta27, beta29, beta3, beta31, beta6, beta8, beta9 
#
# 7 : delta11, delta12, delta13, delta2, delta3, delta4, delta5, delta6, delta7, delta8 
#
# 8 : delta27, delta28, delta34 
#
# 9 : alpha6, alpha7 
#
# 10 : delta16, delta17, delta19




# ------------------------------ New vOTUs / genus-level clusters ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
nrow(t)
# 378


sum(t$n_this_study > 0 & t$n_DBs > 0)
# 102


sum(t$n_DBs == 0)
# 125


G <- unique(t$genus)
G <- G[G != '']
length(G)
# 132


all(t$n_this_study[t$genus == ''] == 0)
# TRUE


N1 <- 0
N2 <- 0
for (x in G) {

    idx <- which(t$genus == x)

    if (all(t$n_DBs[idx] == 0)) { N1 <- N1 + 1 }

    if (any(t$n_this_study[idx] > 0) & any(t$n_DBs[idx] > 0)) { N2 <- N2 + 1 }

}
print(N1)
# 32
print(N2)
# 60




# ------------------------------ Contigs >=10 kb ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/10kb_contigs_number.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

range(t$n_contigs)
# 0 3036




# ------------------------------ Samples with crAss per cohort ------------------------------ #
t <- read.table(
    '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

d <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/NL_crAss_contigs_2.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


for (k in c('LLD', 'LLD2', '300OB', 'IBD')) {

    S <- unique(d$sample_id[d$cohort_name == k])

    mC <- max(sapply(S, function (x) sum(d$sample_id == x)))

    nCs <- length(S)

    nAs <- sum(t$cohort == k)

    cat(k, ':', round(nCs / nAs * 100), '% samples with crAss-like contigs, max', mC, 'crAss-like contigs per sample.\n\n')

}
# LLD : 43 % samples with crAss-like contigs, max 5 crAss-like contigs per sample.
#
# LLD2 : 58 % samples with crAss-like contigs, max 6 crAss-like contigs per sample.
#
# 300OB : 61 % samples with crAss-like contigs, max 6 crAss-like contigs per sample.
#
# IBD : 29 % samples with crAss-like contigs, max 6 crAss-like contigs per sample.




# ------------------------------ crAss in ORFs and 6 frames ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/crAss_contigs_number.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

sum(t$n_crAss)
# 1556

idx <- which((t$n_crAss != t$n_crAss_6frames) | (t$n_crAss != t$n_crAss_orfs))

t[idx, ]
#              sample_id n_crAss n_crAss_6frames n_crAss_orfs
# 263 XXXX       2               1            2
# 1157 XXXX       2               2            1
# 1312 XXXX       1               0            1




# ------------------------------ Properties of crAss contigs ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/NL_crAss_contigs_2.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

round(range(t$length) / 1000, 1)
# 10.0 195.7

round(mean(t$length) / 1000, 1)
# 58.8

sum(t$ends_overlap == 1)
# 255


t <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(t$n_DBs == 0)

length(idx)
# 125

sum(t$translation[idx] == 'TAG|q')
# 49

sum(t$translation[idx] == 'TGA|w')
# 2




# ------------------------------ Reads mapping ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


M <- data.frame(NULL)
for (k in c('LLD', 'LLD2', '300OB', 'IBD')) {

    m <- read.table(
        paste0('/data/umcg-tifn/crAss_analysis/map_reads/', k, '_crAss_abundance.txt'),
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
    )

    if (nrow(M) == 0) { M <- m } else { M <- cbind(M, m) }

}


sele <- t$cl_repres[t$group == '']
all(M[sele, ] == 0)
# TRUE


sele <- t$cl_repres[t$group != '']
V <- apply(M[sele, ], 1, function (v) all(v == 0))
sum(V)
# 17




# ------------------------------ Genome length ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/DATABASES/data_Yutin_2020/TableS1_genomes_used_1.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
range(t$len_nt)
# 8857 192156
range(t$len_nt[t$source == 'cMAG, current study'])
# 61829 192156


d <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/NL_crAss_contigs_2.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
range(d$length)
# 10022 195671
range(d$length[d$ends_overlap == 1])
# 88345 191444


k <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
for (x in c('db673', 'db249', 'db146')) {

    q <- read.table(
        paste0('/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/crAss_in_',x,'/',x,'_ends_overlap.txt'),
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
    )[, 1]

    cat(x, ':', range(k$cl_member_len[k$cl_member %in% q]), '\n')

}
#db673 : 61829 192156
#db249 : 79728 104752
#db146 : 59879 103832




# ------------------------------ 'Complete' vOTU representatives ------------------------------ #
DF <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
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


t <- read.table('/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/NL_crAss_contigs_2.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
V <- t$contig_name[ t$ends_overlap == 1 ]

for (x in c('db146', 'db249', 'db673')) {

    v <- read.table(paste0('/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/crAss_in_', x, '/', x, '_ends_overlap.txt'), sep = '\t', header = FALSE, stringsAsFactors = FALSE)[, 1]
    V <- c(V, v)

}


X <- c()
for (x in DF$cl_repres) {

    y <- k$cl_member[ k$cl_repres == x ]

    if (any(y %in% V)) { X <- c(X, x) }

}
length(X)
# 216


Y <- DF$cl_repres[DF$cl_repres %in% V]

all(Y %in% X)
# TRUE

X[!(X %in% Y)]
# NC_024711_crAssphage

length(Y)
# 215

nrow(DF) - length(Y)
# 163


# genomes crAss001_MH675552 and NC_024711_crAssphage represent vOTUs, and are complete according to the literature, but do not have direct terminal repeats
c('crAss001_MH675552', 'NC_024711_crAssphage') %in% DF$cl_repres
# TRUE TRUE
c('crAss001_MH675552', 'NC_024711_crAssphage') %in% V
# FALSE FALSE




# ------------------------------ Known & novel contigs ------------------------------ #
t <- read.table(
    '/data/umcg-tifn/crAss_analysis/build_tree/CRASS_DB_cl_summary_TRANSL_GROUPS.final.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)


sum(t$n_this_study[(t$n_this_study > 0) & (t$n_DBs > 0)]) / sum(t$n_this_study) * 100
# 75.77121

sum(t$n_this_study[(t$n_this_study > 0) & (t$n_DBs == 0)]) / sum(t$n_this_study) * 100
# 24.22879


t2 <- aggregate(t[, c('n_this_study', 'n_DBs')], by = list(t$genus), FUN = sum)
colnames(t2)[1] <- 'genus'

t2[t2$genus == '', ]
# genus n_this_study n_DBs
#                  0   121

sum(t2$n_this_study[(t2$n_this_study > 0) & (t2$n_DBs > 0)]) / sum(t2$n_this_study) * 100
# 95.50129

sum(t2$n_this_study[(t2$n_this_study > 0) & (t2$n_DBs == 0)]) / sum(t2$n_this_study) * 100
# 4.498715




# ------------------------------ Detected at least once ------------------------------ #
t <- read.table(
    '/home/umcg-agulyaeva/LLD_LLD2_IBD_300OB_assembly/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

M <- read.table(
    '/data/umcg-tifn/crAss_analysis/map_reads/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)
M <- M[grep('\\|g__[^\\|]+$', rownames(M)), ]
rownames(M) <- sub('^.+\\|g__(.+)$', '\\1', rownames(M))


nrow(M)
# 132
sum(apply(M, 1, function (v) any(v > 0)))
# 118


ncol(M)
# 2291
sum(apply(M, 2, function (v) any(v > 0)))
# 1998


for (k in c('LLD', 'LLD2', '300OB', 'IBD')) {

    samples <- t$sample_id[t$cohort == k]
    if (k == 'LLD') { samples <- paste0('X', samples) }

    m <- M[, samples]

    V <- apply(m, 2, function (v) sum(v > 0))
    cat('Mean number of genus-level clusters per ', k, ' sample: ', mean(V), '\n')

}
# Mean number of genus-level clusters per  LLD  sample:  4.25022 
# Mean number of genus-level clusters per  LLD2  sample:  4.721893 
# Mean number of genus-level clusters per  300OB  sample:  5.258389 
# Mean number of genus-level clusters per  IBD  sample:  1.934615




# ------------------------------ longest identified contigs ------------------------------ #
t <- read.table(
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

t$cl_repres <- sapply(t$contig_name, function (x) k$cl_repres[k$cl_member == x])

t[order(t$length, decreasing = T),][1:4, c(4:6, 8)]
# length ends_overlap    contig_name    cl_repres
# 195671            0 NL_crAss001475 OCPX01000081
# 193543            0 NL_crAss001075 OLGJ01000028
# 191444            1 NL_crAss000703 OLGJ01000028
# 190333            1 NL_crAss001463 OGQH01000028




# ------------------------------ trees ------------------------------ #
library(ape)

tr1 <- read.tree('/data/umcg-tifn/crAss_analysis/build_tree/MSA_crAss_TerL_MidpointRooted.2.newick')
Ntip(tr1)
# 371

tr2 <- read.tree('/data/umcg-tifn/crAss_analysis/build_tree/MSA_crAss_portal_MidpointRooted.newick')
Ntip(tr2)
# 371
