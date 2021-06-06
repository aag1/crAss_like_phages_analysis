# based on code by A. Kurilshikov

library(meta)
library(foreach)



ibd  = read.table("IBD_phage_host_abundance_corr.txt", header = T, as.is = T, sep = "\t")
ob   = read.table("300OB_phage_host_abundance_corr.txt", header = T, as.is = T, sep = "\t")
lld  = read.table("LLD_phage_host_abundance_corr.txt", header = T, as.is = T, sep = "\t")
lld2 = read.table("LLD2_phage_host_abundance_corr.txt", header = T, as.is = T, sep = "\t")



ibd_line  = paste(ibd$phage_taxon, ibd$host_taxon)
ob_line   = paste(ob$phage_taxon, ob$host_taxon)
lld_line  = paste(lld$phage_taxon, lld$host_taxon)
lld2_line = paste(lld2$phage_taxon, lld2$host_taxon)

common = unique(c(ibd_line, ob_line, lld_line, lld2_line))

ibd.2  = ibd[match(common, ibd_line), ]
ob.2   = ob[match(common, ob_line), ]
lld.2  = lld[match(common, lld_line), ]
lld2.2 = lld2[match(common, lld2_line), ]

combined = t(data.frame(do.call(cbind, strsplit(common,split = " "))))
combined = data.frame(
                      phage_taxon = combined[,1], host_taxon = combined[,2],
                      lld.R  = lld.2[,3],  lld.P  = lld.2[,4],
                      lld2.R = lld2.2[,3], lld2.P = lld2.2[,4],
                      ibd.R  = ibd.2[,3],  ibd.P  = ibd.2[,4],
                      ob.R   = ob.2[,3],   ob.P   = ob.2[,4]
)



metas = foreach(i = 1:nrow(combined), .combine = rbind) %do% {

  cors = as.vector(combined[i, c(3,7,9)])
  ps = as.vector(combined[i, c(4,8,10)])
  sizes = c(1135,455,298)

  cors2 = cors[!is.na(cors)]
  ps2 = ps[!is.na(cors)]
  sizes2 = sizes[!is.na(cors)]

  if(length(cors2) == 0) {

    c(combined[i,5], combined[i,6])

  } else if (length(cors2) == 1) {

    c(cors2, ps2)

  } else { 

    c1 = metacor(cors2, sizes2, sm = "ZCOR", method.tau = "SJ")

    c(c1$TE.random, c1$pval.random)

  }

}

metas2 = cbind(metas, p.adjust(metas[,2], method = "BH"))

combined = cbind(combined, metas2)

colnames(combined)[11:13] = c("metaR", "metaP", "metaFDR")

combined <- combined[order(combined$metaR, decreasing = T), ]

write.table(
    combined,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    file = 'meta_phage_host_abundance_corr.txt'
)
