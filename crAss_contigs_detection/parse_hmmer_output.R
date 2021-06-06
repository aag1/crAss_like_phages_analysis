.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(optparse)
library(rhmmer)
sessionInfo()



### input parameters
option_list = list(
	make_option('--inF', help='input file (hmmer tblout)'),
	make_option('--outF', help='output file'))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



### read input
t <- read_tblout(opt$inF)
t <- as.data.frame(t, stringsAsFactors = FALSE)

colnames(t)[c(1, 3)] <- c('protein_id', 'profile_id')

t$contig_id <- sub('_[0-9]+$', '', t$protein_id)



### build output table
markers <- c('MCP', 'portal', 'TerL')
contigs <- unique(t$contig_id)

df <- as.data.frame(
        matrix(
            0,
            nrow = length(contigs),
            ncol = length(markers),
            dimnames = list(contigs, markers)),
        stringsAsFactors = FALSE
)


for (x in contigs) {

    for (y in markers) {

        idx <- which(t$contig_id == x & grepl(paste0('^',y,'_'), t$profile_id))

        df[x, y] <- length(idx)

    }

}


markers_num <- apply(df, 1, function (v) sum(v != 0))
sele <- which(markers_num == 3)
df <- df[sele, ]


if (nrow(df) > 0) {

    write.table(
        rownames(df),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        file = opt$outF)

}
