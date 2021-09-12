sessionInfo()




# -------------------- read data --------------------
dat <- data.frame(NULL)

pheno_files <- list.files(path = '../../LLD_LLD2_Info/Pheno_science_imputed_1135/', full.names = T)

for (f in pheno_files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != 'LLDEEPid']

    if (ncol(dat) == 0) { dat <- df } else { dat <- cbind(dat, df) }

}


tab <- read.table(
    'Zhernakova_Table_S1.txt',
    sep = '\t',
    header = T,
    fill = T,
    quote = '',
    comment.char = '',
    stringsAsFactors = F
)




# -------------------- process data --------------------
DF <- data.frame(NULL)

group <- ''

for (i in 1:nrow(tab)) {

    if (tab$Description[i] == '') {

        group <- tolower(tab$Factor.[i])
        if (group == 'smoke') { group <- 'smoking' }
        if (group == 'drugs') { group <- 'medication' }

    } else {

        if (tab$Factor.[i] %in% colnames(dat)) {

            pheno <- tab$Factor.[i]
            transf <- ''

        } else if (paste0(tab$Factor.[i], '_log') %in% colnames(dat)) {

            pheno <- paste0(tab$Factor.[i], '_log')
            transf <- 'log10'

        } else {

            cat('Unrecognized phenotype :', tab$Factor.[i], '\n')

        }


        descr <- tab$Description[i]
        if (descr == 'waist circumference: hip circumference') { descr <- 'waist circumference:hip circumference' }


        units <- tab$Unit[i]
        if (tab$Factor.[i] == 'antrop_gender.F1M2') { units <- '1, female; 2, male' }
        if (units == 'yes') { units <- '0, no; 1, yes' }


        df <- data.frame(
            pheno,
            group,
            descr,
            units,
            transf,
            stringsAsFactors = F
        )

        DF <- rbind(DF, df)

    }

}




# -------------------- write data --------------------
write.table(
    DF,
    quote = F,
    sep = '\t',
    row.names = F,
    file = 'LLD_phenos_list.txt'
)
