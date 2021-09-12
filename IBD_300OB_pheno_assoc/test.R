ids_ibd <- read.table('../../IBD_Info/selected_ibd_samples.txt', stringsAsFactors = FALSE)[, 1]




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

pheno_ibd_2 <- pheno_ibd_2[ids_ibd, ]




idx <- which(pheno_ibd_2$DiagnosisCurrent == 'CD')
df <- pheno_ibd_2[idx, c('DiagnosisCurrent', 'ActiveDisease', 'HarveyBradshaw', 'SSCAI')]
df <- df[order(df$ActiveDisease, -df$HarveyBradshaw), ]
df
cat('\n\n\n')

idx <- which(pheno_ibd_2$DiagnosisCurrent == 'UC')
df <- pheno_ibd_2[idx, c('DiagnosisCurrent', 'ActiveDisease', 'HarveyBradshaw', 'SSCAI')]
df <- df[order(df$ActiveDisease, -df$SSCAI), ]
df
cat('\n\n\n')

idx <- which(!(pheno_ibd_2$DiagnosisCurrent %in% c('CD', 'UC')))
df <- pheno_ibd_2[idx, c('DiagnosisCurrent', 'ActiveDisease', 'HarveyBradshaw', 'SSCAI')]
df <- df[order(df$ActiveDisease, -df$SSCAI), ]
df
cat('\n\n\n')

any(is.na(pheno_ibd_2$DiagnosisCurrent))
