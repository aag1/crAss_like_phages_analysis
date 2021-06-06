.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(ape)
library(phangorn)
sessionInfo()



################################################## ASSIGN GROUPS ##################################################
### read and root tree
tr <- read.tree('MSA_crAss_TerL.fasta.treefile')

tr <- midpoint(tr)
tr <- read.tree(text = write.tree(tr))

write.tree(tr, file = 'MSA_crAss_TerL_MidpointRooted.1.newick')



### read data
DF <- read.table(
    'CRASS_DB_cl_summary_TRANSL.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(DF) <- DF$cl_repres


t <- read.table(
    '/data/umcg-tifn/crAss_analysis/crAss_contigs_detection/CRASS_DB_cl.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
rownames(t) <- t$cl_member


d <- read.table(
    '/data/umcg-tifn/DATABASES/data_Yutin_2020/TableS1_genomes_used_1.txt',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
) 



### assign groups and subgroups
# 1
t$group <- ''
t[rownames(d), 'group'] <- d$group

t$subgroup <- ''
t[rownames(d), 'subgroup'] <- d$subgroup


# 2
DF$group <- ''
DF$subgroup <- ''

for (x in DF$cl_repres) {

    G <- t$group[ t$cl_repres == x ]

    G <- G[ G != '' ]
    if (length(G) == 0) { next }

    G <- unique(G)
    if (length(G) == 1) { DF[x, 'group'] <- G } else {

        DF[x, 'group'] <- t$group[ t$cl_member == x ]
        warning('vOTU represented by ', x, ' contains sequences belonging to multiple groups: ', paste(G, collapse = ','), '!')

    }


    S <- t$subgroup[ t$cl_repres == x ]

    S <- S[ S != '' ]
    if (length(S) == 0) { next }

    S <- unique(S)
    if (length(S) == 1) { DF[x, 'subgroup'] <- S } else {

        DF[x, 'subgroup'] <- t$subgroup[ t$cl_member == x ]
        warning('vOTU represented by ', x, ' contains sequences belonging to multiple subgroups: ', paste(S, collapse = ','), '!')

    }

}


# 3
G <- c('alpha_gamma', 'beta', 'delta', 'epsilon', 'zeta')

for (g in G) {

    X <- rownames(DF)[ DF$group == g ]
    X <- X[X %in% tr$tip.label]

    mrca <- getMRCA(tr, tip = X)

    mrca_desc <- Descendants(tr, node = mrca, type = 'tips')[[1]]
    mrca_desc <- tr$tip.label[mrca_desc]

    if (all(DF[mrca_desc, 'group'] %in% c(g, ''))) {

        DF[mrca_desc, 'group'] <- g

    } else {

        z <- unique(DF[mrca_desc, 'group'])
        z <- z[z != '']
        warning('There are representatives of multile groups among the ', g, ' group MRCA descendants: ', paste(z, collapse = ','), '!')

    }

}


S <- unique(DF$subgroup)
S <- S[ S != '' ]

for (s in S) {

    X <- rownames(DF)[ DF$subgroup == s ]
    X <- X[X %in% tr$tip.label]
    if (length(X) == 1) { next }

    mrca <- getMRCA(tr, tip = X)

    mrca_desc <- Descendants(tr, node = mrca, type = 'tips')[[1]]
    mrca_desc <- tr$tip.label[mrca_desc]

    if (all(DF[mrca_desc, 'subgroup'] %in% c(s, ''))) {

        DF[mrca_desc, 'subgroup'] <- s

    } else {

        z <- unique(DF[mrca_desc, 'subgroup'])
        z <- z[z != '']
        warning('There are representatives of multile subgroups among the ', s, ' subgroup MRCA descendants: ', paste(z, collapse = ','), '!')

    }
}


# 4
DF$group[DF$group %in% c('env', 'env_beta', 'env_delta', 'env_epsilon', 'outgroup', 'Flavob_phages')] <- ''
DF$subgroup[DF$subgroup %in% c('Fpv3', 'phi14_2')] <- ''



################################################## ASSIGN VCs, GENERA ##################################################
source('/home/umcg-agulyaeva/crAss_analysis/crAss_contigs_detection/function_import_cl_results.R')



### read data
tab1 <- read.table(
    '/data/umcg-tifn/crAss_analysis/run_vcontact/genome_by_genome_overview.csv',
    sep = ',',
    header = TRUE,
    row.names = 2,
    stringsAsFactors = FALSE
)


tab2 <- import_cl_results('/data/umcg-tifn/crAss_analysis/genus_level_clusters/CRASS_DB_cl_0-50.clstr')
rownames(tab2) <- tab2$cl_member



### add vConTACT2 data
DF$VC.Status <- tab1[rownames(DF), 'VC.Status']
DF$VC.Subcluster <- tab1[rownames(DF), 'VC.Subcluster']



### add clusters at 50% nucleotide similarity
DF$genus_raw <- tab2[DF$cl_repres, 'cl_repres']



### write table
write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'CRASS_DB_cl_summary_TRANSL_GROUPS.1.txt'
)
