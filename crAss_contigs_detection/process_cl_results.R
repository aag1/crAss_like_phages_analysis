sessionInfo()
source('/home/umcg-agulyaeva/crAss_analysis/crAss_contigs_detection/function_import_cl_results.R')



### complete genomes from different databases
L <- list()

tab <- read.table('NL_crAss_contigs.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
L[['CRASS_DB']] <- tab$contig_name[ tab$ends_overlap == 1 ]

for (x in c('db146', 'db249', 'db673')) {

    L[[x]] <- read.table(paste0('crAss_in_', x, '/', x, '_ends_overlap.txt'), sep = '\t', header = FALSE, stringsAsFactors = FALSE)[, 1]

}



### raw clustering results
DF <- import_cl_results('CRASS_DB_95-85.clstr')



### select cluster representative
for (rep in unique(DF$cl_repres)) {

    idx <- which(DF$cl_repres == rep)

    memb <- DF$cl_member[ idx ]

    new_rep <- rep


    for (x in names(L)) {

        sele <- which(memb %in% L[[x]])
        
        if (length(sele) > 0)  { new_rep <- memb[ sele ][1] }

    }


    for (x in c('IAS_virus_KJ003983', 'crAss001_MH675552', 'NC_024711_crAssphage')) {

        if (x %in% memb) { new_rep <- x }

    }


    DF$cl_repres[ idx ] <- new_rep

}



### write table
write.table(
    DF,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'CRASS_DB_cl.txt'
)
