# based on code by R.A.A.A. Ruigrok

.libPaths('/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB')
library(foreach)
sessionInfo()




# -------------------- read data --------------------
# samples per cohort
k <- read.table('../from_Peregrine/sample_cohort.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

ids_lld_0 <- k$sample_id[k$cohort == 'LLD']
ids_lld   <- paste0('X', ids_lld_0)

ids_300ob <- k$sample_id[k$cohort == '300OB']

ids_ibd <- read.table('../../IBD_Info/selected_ibd_samples.txt', stringsAsFactors = FALSE)[, 1]



# crAss abundance
crass <- read.table(
    '../from_Peregrine/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

crass <- crass[grep('\\|s__', rownames(crass), invert = TRUE), ]  # no vOTUs

rownames(crass) <- sub('^o__', '', rownames(crass))  # no order/family/genus designations
rownames(crass) <- gsub('\\|[a-z]__', '|', rownames(crass))

crass_lld <- crass[, ids_lld]
sele <- rownames(crass_lld)[ rowSums(crass_lld > 0) > length(ids_lld) * 0.05 ]  # select taxanomic clusters present in >5% LLD samples
crass <- crass[sele, ]

crass <- t(crass)



# microbial abundance
metaphlan = read.table(
    '../from_Peregrine/LLD_LLD2_300OB_IBD_merged_abundance_table.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

metaphlan <- metaphlan[, colnames(metaphlan) != 'NCBI_tax_id']

colnames(metaphlan) <- sub('_metaphlan$', '', colnames(metaphlan))

metaphlan <- t(metaphlan)



# LLD phenotypes
pheno_lld <- data.frame(NULL)

pheno_files <- list.files(path = '../../LLD_LLD2_Info/Pheno_science_imputed_1135/', full.names = TRUE)

for (f in pheno_files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != 'LLDEEPid']

    if (ncol(pheno_lld) == 0) { pheno_lld <- df } else { pheno_lld <- cbind(pheno_lld, df) }

}

key_lld <- read.table('../../LLD_LLD2_Info/LLD_GTMicrob.txt', sep = '\t', stringsAsFactors = FALSE)
key_lld$V2 <- paste0('X', key_lld$V2)
rownames(pheno_lld) <- sapply(rownames(pheno_lld), function (x) key_lld$V2[ key_lld$V1 == x ])

sum(ids_lld %in% rownames(pheno_lld)) # 1135

any(is.na(pheno_lld[, 'antrop_age']))         # FALSE
any(is.na(pheno_lld[, 'antrop_gender.F1M2'])) # FALSE



# IBD phenotypes
pheno_ibd_1 <- read.table(
    '../../IBD_Info/data_from_Renate/LLDIBDMIBS_biomarkers_raw.csv',
    sep = ',',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(pheno_ibd_1) %in% ids_lld_0)
rownames(pheno_ibd_1)[idx] <- paste0('X', rownames(pheno_ibd_1)[idx])

sum(ids_lld %in% rownames(pheno_ibd_1)) # 1129
sum(ids_ibd %in% rownames(pheno_ibd_1)) # 431


pheno_ibd_2 <- read.table(
    '../../IBD_Info/data_from_Renate/LLD_IBD_meta_201020.txt',
    quote = '',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

key_ibd <- read.table('../../IBD_Info/rename_IBD.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(pheno_ibd_2) <- sapply(rownames(pheno_ibd_2), function (x) ifelse(x %in% key_ibd$Classic[key_ibd$Classic != 'XXXX'], key_ibd$old[key_ibd$Classic == x], x))
rownames(pheno_ibd_2)[rownames(pheno_ibd_2) == 'XXXX'] <- 'YYYY'

idx <- which(rownames(pheno_ibd_2) %in% key_lld$V1)
rownames(pheno_ibd_2)[idx] <- sapply(rownames(pheno_ibd_2)[idx], function (x) key_lld$V2[ key_lld$V1 == x ])

sum(ids_lld %in% rownames(pheno_ibd_2)) # 1135
sum(ids_ibd %in% rownames(pheno_ibd_2)) # 455

any(is.na(pheno_ibd_2[, 'AgeAtFecalSampling'])) # FALSE
any(is.na(pheno_ibd_2[, 'Sex']))                # FALSE



# 300OB phenotypes
key_300ob <- read.table('../../300OB_Info/key_300OB.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


pheno_300ob_1 <- read.table(
    '../../300OB_Info/300OB_65phenotype.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(pheno_300ob_1) %in% key_300ob$ID)
rownames(pheno_300ob_1)[idx] <- sapply(rownames(pheno_300ob_1)[idx], function (x) key_300ob$G_id[ key_300ob$ID == x ])

sum(ids_300ob %in% rownames(pheno_300ob_1)) # 298

mean(pheno_300ob_1$Length[pheno_300ob_1$sex == 1])  # 164.4259
mean(pheno_300ob_1$Length[pheno_300ob_1$sex == 2])  # 177.2186
mean(pheno_300ob_1$Hip_circumference[pheno_300ob_1$sex == 1])  # 113.4222
mean(pheno_300ob_1$Hip_circumference[pheno_300ob_1$sex == 2])  # 109.9281

any(is.na(pheno_300ob_1[, 'age'])) # FALSE
any(is.na(pheno_300ob_1[, 'sex'])) # FALSE


pheno_300ob_2 <- read.table(
    '../../300OB_Info/4tb8spry4b-1/300OB_metabolicSyndrome.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(pheno_300ob_2) %in% key_300ob$ID)
rownames(pheno_300ob_2)[idx] <- sapply(rownames(pheno_300ob_2)[idx], function (x) key_300ob$G_id[ key_300ob$ID == x ])

sum(ids_300ob %in% rownames(pheno_300ob_2)) # 298




# -------------------- LLD vs. IBD --------------------
sele <- grep('\\|g__(Bacteroides|Prevotella|Porphyromonas|Parabacteroides)$', colnames(metaphlan), value = T)

sc <- min(metaphlan[metaphlan > 0]) / 2     # small constant

hosts <- as.data.frame(log10(metaphlan[, sele] + sc))    # log transformation to get normally distributed data



IDS <- c(ids_lld, ids_ibd)

Cohort <- c(rep(0, length(ids_lld)), rep(1, length(ids_ibd)))

names(Cohort) <- IDS



DF1 <- cbind(
    data.frame(Cohort),
    pheno_ibd_2[IDS, c('AgeAtFecalSampling', 'Sex', 'SUMOFKHTOT', 'group_meat', 'how_often_coffee', 'laxatives')],
    pheno_ibd_1[IDS, c('ChrA', 'HBD2')]
)

DF1[sapply(DF1, is.character)] <- lapply(DF1[sapply(DF1, is.character)], as.factor)

DF2 <- cbind(DF1, hosts[IDS, ])



ids_lld_sele <- ids_lld[ apply(DF1[ids_lld, ], 1, function (v) all(!is.na(v))) ]
length(ids_lld_sele) # 1098

ids_ibd_sele <- ids_ibd[ apply(DF1[ids_ibd, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_sele) # 409



OUT1 <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF1, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        LLD_positive_samples = round(sum(crass[ids_lld, i] > 0) / length(ids_lld) * 100, 2),
        IBD_positive_samples = round(sum(crass[ids_ibd, i] > 0) / length(ids_ibd) * 100, 2),
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}



OUT2 <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF2, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}



OUT = data.frame(
    crAss_taxon = OUT1$crAss_taxon,

    LLD_positive_samples = OUT1$LLD_positive_samples,
    IBD_positive_samples = OUT1$IBD_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'fdr'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'fdr')
)

OUT = OUT[order(OUT$FDR.glm1), ]

write.table(
    OUT,
    file = 'LLD_vs_IBD.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- IBD: CD vs. UC --------------------
Diagnosis <- setNames(pheno_ibd_2$DiagnosisCurrent, rownames(pheno_ibd_2))

Diagnosis <- Diagnosis[ids_ibd]

Diagnosis <- Diagnosis[!is.na(Diagnosis)]

unique(Diagnosis) # CD UC IBDU ReconsideringDiagnosis IBDI MicroscopicColitis

ids_ibd_cd <- names(Diagnosis)[Diagnosis == 'CD']
ids_ibd_uc <- names(Diagnosis)[Diagnosis == 'UC']

length(ids_ibd_cd) # 256
length(ids_ibd_uc) # 171

IDS <- c(ids_ibd_cd, ids_ibd_uc)

Diagnosis <- Diagnosis[IDS]

Diagnosis <- ifelse(Diagnosis == 'CD', 0, 1)



DF <- cbind(
    data.frame(Diagnosis),
    pheno_ibd_2[IDS, c('AgeAtFecalSampling', 'Sex', 'SUMOFKHTOT', 'group_meat', 'how_often_coffee', 'laxatives')],
    pheno_ibd_1[IDS, c('ChrA', 'HBD2')]
)

DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)



ids_ibd_cd_sele <- ids_ibd_cd[ apply(DF[ids_ibd_cd, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_cd_sele) # 229

ids_ibd_uc_sele <- ids_ibd_uc[ apply(DF[ids_ibd_uc, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_uc_sele) # 154



OUT <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        CD_positive_samples = round(sum(crass[ids_ibd_cd, i] > 0) / length(ids_ibd_cd) * 100, 2),
        UC_positive_samples = round(sum(crass[ids_ibd_uc, i] > 0) / length(ids_ibd_uc) * 100, 2),
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}

OUT$FDR <- p.adjust(OUT$P, method = 'fdr')

OUT = OUT[order(OUT$FDR), ]

write.table(
    OUT,
    file = 'IBD_CD_vs_UC.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- IBD: exclusively colonic vs. ileal-inclusive disease location --------------------
Location <- setNames(pheno_ibd_2$DiseaseLocation, rownames(pheno_ibd_2))

Location <- Location[ids_ibd]

Location <- Location[!is.na(Location)]

unique(Location) # ileum colon both

ids_ibd_ec <- names(Location)[Location == 'colon']
ids_ibd_ii <- names(Location)[Location %in% c('ileum', 'both')]

length(ids_ibd_ec) # 213
length(ids_ibd_ii) # 195

IDS <- c(ids_ibd_ec, ids_ibd_ii)

Location <- Location[IDS]

Location <- ifelse(Location == 'colon', 0, 1)



DF <- cbind(
    data.frame(Location),
    pheno_ibd_2[IDS, c('AgeAtFecalSampling', 'Sex', 'SUMOFKHTOT', 'group_meat', 'how_often_coffee', 'laxatives')],
    pheno_ibd_1[IDS, c('ChrA', 'HBD2')]
)

DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)



ids_ibd_ec_sele <- ids_ibd_ec[ apply(DF[ids_ibd_ec, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_ec_sele) # 197

ids_ibd_ii_sele <- ids_ibd_ii[ apply(DF[ids_ibd_ii, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_ii_sele) # 170



OUT <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        ec_positive_samples = round(sum(crass[ids_ibd_ec, i] > 0) / length(ids_ibd_ec) * 100, 2),
        ii_positive_samples = round(sum(crass[ids_ibd_ii, i] > 0) / length(ids_ibd_ii) * 100, 2),
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}

OUT$FDR <- p.adjust(OUT$P, method = 'fdr')

OUT = OUT[order(OUT$FDR), ]

write.table(
    OUT,
    file = 'IBD_excl_colonic_vs_ileum_incl.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- IBD: remission vs. active disease --------------------
Activity <- setNames(pheno_ibd_2$ActiveDisease, rownames(pheno_ibd_2))

Activity <- Activity[ids_ibd]

Activity <- Activity[!is.na(Activity)]

unique(Activity) # NotActive Active

ids_ibd_re <- names(Activity)[Activity == 'NotActive']
ids_ibd_ac <- names(Activity)[Activity == 'Active']

length(ids_ibd_re) # 338
length(ids_ibd_ac) # 110

IDS <- c(ids_ibd_re, ids_ibd_ac)

Activity <- Activity[IDS]

Activity <- ifelse(Activity == 'NotActive', 0, 1)



DF <- cbind(
    data.frame(Activity),
    pheno_ibd_2[IDS, c('AgeAtFecalSampling', 'Sex', 'SUMOFKHTOT', 'group_meat', 'how_often_coffee', 'laxatives')],
    pheno_ibd_1[IDS, c('ChrA', 'HBD2')]
)

DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)



ids_ibd_re_sele <- ids_ibd_re[ apply(DF[ids_ibd_re, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_re_sele) # 301

ids_ibd_ac_sele <- ids_ibd_ac[ apply(DF[ids_ibd_ac, ], 1, function (v) all(!is.na(v))) ]
length(ids_ibd_ac_sele) # 102



OUT <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        re_positive_samples = round(sum(crass[ids_ibd_re, i] > 0) / length(ids_ibd_re) * 100, 2),
        ac_positive_samples = round(sum(crass[ids_ibd_ac, i] > 0) / length(ids_ibd_ac) * 100, 2),
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}

OUT$FDR <- p.adjust(OUT$P, method = 'fdr')

OUT = OUT[order(OUT$FDR), ]

write.table(
    OUT,
    file = 'IBD_remission_vs_active_disease.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- LLD vs. 300OB --------------------
IDS <- c(ids_lld, ids_300ob)

DF1 <- data.frame(
    Cohort = c(
        rep(0, length(ids_lld)),
        rep(1, length(ids_300ob))),
    Age = c(
        pheno_lld[ids_lld, 'antrop_age'],
        pheno_300ob_1[ids_300ob, 'age']),
    Gender = c(
        pheno_lld[ids_lld, 'antrop_gender.F1M2'],
        pheno_300ob_1[ids_300ob, 'sex'])
)

rownames(DF1) <- IDS

DF1[sapply(DF1, is.character)] <- lapply(DF1[sapply(DF1, is.character)], as.factor)

DF2 <- cbind(DF1, hosts[IDS, ])



OUT1 <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF1, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        LLD_positive_samples = round(sum(crass[ids_lld, i] > 0) / length(ids_lld) * 100, 2),
        OB_positive_samples = round(sum(crass[ids_300ob, i] > 0) / length(ids_300ob) * 100, 2),
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}



OUT2 <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[IDS, i] > 0 ~ ., data = DF2, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}



OUT = data.frame(
    crAss_taxon = OUT1$crAss_taxon,

    LLD_positive_samples = OUT1$LLD_positive_samples,
    OB_positive_samples = OUT1$OB_positive_samples,

    Beta.glm1 = OUT1$Beta,
    SE.glm1   = OUT1$SE,
    Z.glm1    = OUT1$Z,
    P.glm1    = OUT1$P,
    FDR.glm1  = p.adjust(OUT1$P, method = 'fdr'),

    Beta.glm2 = OUT2$Beta,
    SE.glm2   = OUT2$SE,
    Z.glm2    = OUT2$Z,
    P.glm2    = OUT2$P,
    FDR.glm2  = p.adjust(OUT2$P, method = 'fdr')
)

OUT = OUT[order(OUT$FDR.glm1), ]

write.table(
    OUT,
    file = 'LLD_vs_300OB.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)




# -------------------- 300OB: with vs. without metabolic syndrome --------------------
Syndrome <- setNames(pheno_300ob_2$MetS.NCEPsom, rownames(pheno_300ob_2))

Syndrome <- Syndrome[ids_300ob]

any(is.na(Syndrome)) # FALSE

unique(Syndrome) # 3 5 2 4 1 0

Syndrome <- ifelse(Syndrome >= 3, 1, 0)

ids_300ob_no <- names(Syndrome)[Syndrome == 0]
ids_300ob_sy <- names(Syndrome)[Syndrome == 1]

length(ids_300ob_no) # 137
length(ids_300ob_sy) # 161



DF <- cbind(
    data.frame(Syndrome),
    pheno_300ob_1[ids_300ob, c('age', 'sex')]
)

DF[sapply(DF, is.character)] <- lapply(DF[sapply(DF, is.character)], as.factor)



OUT <- foreach (i = 1:ncol(crass), .combine = rbind) %do% {

    glm.sum <- summary(glm(crass[ids_300ob, i] > 0 ~ ., data = DF, family = binomial(link = 'logit')))

    data.frame(
        crAss_taxon = colnames(crass)[i],
        no_positive_samples = round(sum(crass[ids_300ob_no, i] > 0) / length(ids_300ob_no) * 100, 2),
        sy_positive_samples = round(sum(crass[ids_300ob_sy, i] > 0) / length(ids_300ob_sy) * 100, 2),
        Beta = glm.sum$coef[2,1],
        SE = glm.sum$coef[2,2],
        Z = glm.sum$coef[2,3],
        P = glm.sum$coef[2,4]
    )

}

OUT$FDR <- p.adjust(OUT$P, method = 'fdr')

OUT = OUT[order(OUT$FDR), ]

write.table(
    OUT,
    file = '300OB_absence_vs_presence_metabolic_syndrome.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)
