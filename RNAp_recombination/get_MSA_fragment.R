.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(bio3d)
sessionInfo()



MSA <- read.fasta('delta27_genomes_circular_MSA.fasta')$ali
idx <- which(MSA['OLNE01000081', ] != '-')



# MSA fragment corresponding to the OLNE01000081 RNAP and VP02740
# see /data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/crAss_key_doms_in_OLNE01000081_TAGq.txt
# see /data/umcg-tifn/crAss_analysis/genome_maps/ANNOTATIONS/OLNE01000081/OLNE01000081_TAGq_AA_coordinates.txt
coo1 <- 50656   # OLNE01000081_39 5'-terminus
coo2 <- 72821   # OLNE01000081_40 3'-terminus

coo1 <- idx[coo1]
coo2 <- idx[coo2]



write.fasta(
    ids = paste0(rownames(MSA), '_RNAP_VP02740'),
    seqs = MSA[, coo1:coo2],
    file = 'delta27_RNAP_VP02740_MSA.fasta'
)
