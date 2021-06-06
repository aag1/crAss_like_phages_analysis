.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(rhmmer)
library(IRanges)
sessionInfo()



### input parameters
option_list = list(
	make_option('--domF', help='list of domains'),
	make_option('--inF', help='input file (hmmer domtblout)'),
	make_option('--outF', help='output file'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read input
key_domains <- read.table(opt$domF, sep = '\t', header = FALSE, stringsAsFactors = FALSE)[, 1]

t <- read_domtblout(opt$inF)
t <- as.data.frame(t, stringsAsFactors = FALSE)

colnames(t)[c(1, 4)] <- c('protein_id', 'profile_id')
t$key_domain <- unlist(lapply(strsplit(t$profile_id, '_'), function (v) v[1]))



### build output table
proteins <- unique(t$protein_id)

df <- as.data.frame(
        matrix(
            '-',
            nrow = length(proteins),
            ncol = length(key_domains),
            dimnames = list(proteins, key_domains)
        ),
        stringsAsFactors = FALSE
)

for (p in proteins) {
    for (k in key_domains) {

        idx <- which(t$protein_id == p & grepl(paste0('^',k,'_'), t$profile_id))
        if (length(idx) == 0) { next }

        ir <- IRanges(start=t$env_from[idx], end=t$env_to[idx])
        cl <- reduce(ir)
        cl <- as.data.frame(cl)
        cl$coo <- paste0(cl$start, '-', cl$end)

        df[p, k] <- paste(cl$coo, collapse = ';')

    }
}

write.table(df, sep = '\t', quote = FALSE, file = opt$outF)
