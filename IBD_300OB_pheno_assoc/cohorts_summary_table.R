sessionInfo()




DF <- data.frame(matrix(NA, ncol = 6, nrow = 4), stringsAsFactors = FALSE)
rownames(DF) <- c('LLD', 'LLD2', '300OB', 'IBD')
colnames(DF) <- c('n_samples', 'M_pct', 'F_pct', 'age', 'age_range', 'BMI')

t <- read.table(
    '../from_Peregrine/sample_cohort.txt',
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)

NN <- function (x) { format(round(x, 1), nsmall = 1) }




### LLD
ids <- t$sample_id[t$cohort == 'LLD']
ids <- paste0('X', ids)



PH <- data.frame(NULL)

files <- list.files(path = '../../LLD_LLD2_Info/Pheno_science_imputed_1135/', full.names = TRUE)

for (f in files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != 'LLDEEPid']

    if (ncol(PH) == 0) { PH <- df } else { PH <- cbind(PH, df) }

}

k <- read.table('../../LLD_LLD2_Info/LLD_GTMicrob.txt', sep = '\t', stringsAsFactors = FALSE)
k$V2 <- paste0('X', k$V2)
rownames(PH) <- sapply(rownames(PH), function (x) k$V2[ k$V1 == x ])

all(ids %in% rownames(PH))
PH <- PH[ids, ]



DF['LLD', 'n_samples'] <- length(ids)

DF['LLD', 'M_pct'] <- NN(sum(PH[, 'antrop_gender.F1M2'] == 2) / DF['LLD', 'n_samples'] * 100)
DF['LLD', 'F_pct'] <- NN(sum(PH[, 'antrop_gender.F1M2'] == 1) / DF['LLD', 'n_samples'] * 100)

DF['LLD', 'age'] <- paste(NN(mean(PH[, 'antrop_age'])), '±', NN(sd(PH[, 'antrop_age'])))
DF['LLD', 'age_range'] <- paste(NN(min(PH[, 'antrop_age'])), '-', NN(max(PH[, 'antrop_age'])))

DF['LLD', 'BMI'] <- paste(NN(mean(PH[, 'antrop_BMI'])), '±', NN(sd(PH[, 'antrop_BMI'])))




### LLD2
ids <- t$sample_id[t$cohort == 'LLD2']



k <- read.table(
    '../../LLD_LLD2_Info/key_LLD_baseline_fup_338sample_pairs.txt',
    sep = '\t',
    row.names = 3,
    header = TRUE,
    stringsAsFactors = FALSE
)

PH <- read.table(
    '../../LLD_LLD2_Info/LLD2_Phenotypes/data_pheno_LLD_base_fup_338pairs_62pheno_18med_17disease_log_min_10participants.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

PH <- PH[ grepl('_F$', rownames(PH)), ]

rownames(PH) <- sapply(rownames(PH), function (x) k$barcode_fup[k$lld_id_fup == x])
rownames(PH) <- sub('\\.bam$', '', rownames(PH))

all(ids %in% rownames(PH))
PH <- PH[ids, ]



DF['LLD2', 'n_samples'] <- length(ids)

DF['LLD2', 'M_pct'] <- NN(sum(PH[, 'antrop_gender.F1M2'] == 2) / DF['LLD2', 'n_samples'] * 100)
DF['LLD2', 'F_pct'] <- NN(sum(PH[, 'antrop_gender.F1M2'] == 1) / DF['LLD2', 'n_samples'] * 100)

DF['LLD2', 'age'] <- paste(NN(mean(PH[, 'antrop_age'])), '±', NN(sd(PH[, 'antrop_age'])))
DF['LLD2', 'age_range'] <- paste(NN(min(PH[, 'antrop_age'])), '-', NN(max(PH[, 'antrop_age'])))

DF['LLD2', 'BMI'] <- paste(NN(mean(PH[, 'antrop_BMI'])), '±', NN(sd(PH[, 'antrop_BMI'])))




### 300OB
ids <- t$sample_id[t$cohort == '300OB']



k <- read.table('../../300OB_Info/key_300OB.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

PH <- read.table(
    '../../300OB_Info/300OB_65phenotype.txt',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

idx <- which(rownames(PH) %in% k$ID)
rownames(PH)[idx] <- sapply(rownames(PH)[idx], function (x) k$G_id[ k$ID == x ])

all(ids %in% rownames(PH))
PH <- PH[ids, ]



DF['300OB', 'n_samples'] <- length(ids)

DF['300OB', 'M_pct'] <- NN(sum(PH[, 'sex'] == 2) / DF['300OB', 'n_samples'] * 100)
DF['300OB', 'F_pct'] <- NN(sum(PH[, 'sex'] == 1) / DF['300OB', 'n_samples'] * 100)

DF['300OB', 'age'] <- paste(NN(mean(PH[, 'age'])), '±', NN(sd(PH[, 'age'])))
DF['300OB', 'age_range'] <- paste(NN(min(PH[, 'age'])), '-', NN(max(PH[, 'age'])))

DF['300OB', 'BMI'] <- paste(NN(mean(PH[, 'BMI'])), '±', NN(sd(PH[, 'BMI'])))




### IBD
ids <- read.table('../../IBD_Info/selected_ibd_samples.txt', stringsAsFactors = FALSE)[, 1]



PH <- read.table(
    '../../IBD_Info/data_from_Renate/LLD_IBD_meta_201020.txt',
    quote = '',
    sep = '\t',
    row.names = 1,
    header = TRUE,
    stringsAsFactors = FALSE
)

k <- read.table('../../IBD_Info/rename_IBD.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(PH) <- sapply(rownames(PH), function (x) ifelse(x %in% k$Classic[k$Classic != 'XXXX'], k$old[k$Classic == x], x))
rownames(PH)[rownames(PH) == 'XXXX'] <- 'YYYY'

all(ids %in% rownames(PH))
PH <- PH[ids, ]



DF['IBD', 'n_samples'] <- length(ids)

DF['IBD', 'M_pct'] <- NN(sum(PH[, 'Sex'] == 'male') / DF['IBD', 'n_samples'] * 100)
DF['IBD', 'F_pct'] <- NN(sum(PH[, 'Sex'] == 'female') / DF['IBD', 'n_samples'] * 100)

DF['IBD', 'age'] <- paste(NN(mean(PH[, 'AgeAtFecalSampling'])), '±', NN(sd(PH[, 'AgeAtFecalSampling'])))
DF['IBD', 'age_range'] <- paste(NN(min(PH[, 'AgeAtFecalSampling'])), '-', NN(max(PH[, 'AgeAtFecalSampling'])))

sum(is.na(PH[, 'BMI']))
DF['IBD', 'BMI'] <- paste(NN(mean(PH[, 'BMI'], na.rm = TRUE)), '±', NN(sd(PH[, 'BMI'], na.rm = TRUE)))




write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'cohorts_summary_table.txt'
)
