# based on code by A. Kurilshikov

.libPaths("/groups/umcg-tifn/tmp01/users/umcg-agulyaeva/SOFTWARE/R_LIB")
library(foreach)
sessionInfo()




# -------------------- read data --------------------
crass = read.table("../from_Peregrine/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style.txt", header = T, as.is = T)



linkage_file = read.table("../../LLD_LLD2_Info/LLD_GTMicrob.txt", as.is = T)

linkage_file$V2 = paste0("X", linkage_file$V2)

shared_samples = intersect(linkage_file$V2, colnames(crass))

linkage_file2 = linkage_file[match(shared_samples, linkage_file$V2), ]
# 'match' returns a vector of the positions of (first) matches of its first argument in its second



crass1 = crass[grep("\\|s__", rownames(crass), invert = T), ]

crass2 = t(crass1[, match(shared_samples, colnames(crass1))])

crass3 = crass2[, colSums(crass2 > 0) > nrow(crass2) * 0.05]



pheno <- data.frame(NULL)

pheno_files <- list.files(path = "../../LLD_LLD2_Info/Pheno_science_imputed_1135/", full.names = T)

for (f in pheno_files) {

    df <- read.table(f, header = T, as.is = T)
    df <- df[, colnames(df) != "LLDEEPid"]

    if (ncol(pheno) == 0) { pheno <- df } else { pheno <- cbind(pheno, df) }

}

pheno2 = pheno[match(linkage_file2$V1, rownames(pheno)), ]



metaphlan = read.table("../from_Peregrine/LLD_LLD2_300OB_IBD_merged_abundance_table.txt", sep = "\t", row.names = 1, header = T, as.is = T)

metaphlan = metaphlan[, colnames(metaphlan) != "NCBI_tax_id"]

colnames(metaphlan) = sub("_metaphlan$", "", colnames(metaphlan))

metaphlan = t(metaphlan)

metaphlan = metaphlan[rownames(crass3), ]




# -------------------- Spearman correlation: crAss abundance --------------------
result.quant = foreach(i = 1:ncol(crass3), .combine = rbind) %:%

  foreach(j = 1:ncol(pheno2), .combine = rbind) %do% {

    c1 = cor.test(pheno2[, j], crass3[, i], method = "spearman")

    data.frame(
        pheno = colnames(pheno2)[j],
        taxon = colnames(crass3)[i],
        Rsp = c1$estimate,
        P = c1$p.value
    )

  }

warnings()




# -------------------- Spearman correlation: crAss presence/absence --------------------
result.bin = foreach(i = 1:ncol(crass3), .combine = rbind) %:%

  foreach(j = 1:ncol(pheno2), .combine = rbind) %do% {

    c1 = cor.test(pheno2[, j], as.integer(crass3[, i] > 0), method = "spearman")

    data.frame(
        pheno = colnames(pheno2)[j],
        taxon = colnames(crass3)[i],
        Rsp = c1$estimate,
        P = c1$p.value
    )

  }

warnings()




# -------------------- Logistic regression: crAss presence/absence --------------------
result.glm1 = foreach(i = 1:ncol(crass3), .combine = rbind) %:%

  foreach(j = 1:ncol(pheno2), .combine = rbind) %do% {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, j], family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )
    
  }

warnings()




# -------------------- Logistic regression: crAss presence/absence; adjusted for Age & Sex --------------------
result.glm2 = foreach(i = 1:ncol(crass3), .combine = rbind) %:%

  foreach(j = 1:ncol(pheno2), .combine = rbind) %do% {

    if (colnames(pheno2)[j] == "antrop_age") {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, "antrop_age"] + pheno2[, "antrop_gender.F1M2"], family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    } else if (colnames(pheno2)[j] == "antrop_gender.F1M2") {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, "antrop_age"] + pheno2[, "antrop_gender.F1M2"], family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[3,1],
            SE = c1$coef[3,2],
            Z = c1$coef[3,3],
            P = c1$coef[3,4]
        )

    } else {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, j] + pheno2[, "antrop_age"] + pheno2[, "antrop_gender.F1M2"], family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    }
    
  }

warnings()




# -------------------- Logistic regression: crAss presence/absence; adjusted for Age, Sex & Hosts --------------------
sele <- grep("\\|g__(Bacteroides|Prevotella|Porphyromonas|Parabacteroides)$", colnames(metaphlan), value = T)

sc <- min(metaphlan[metaphlan > 0]) / 2     # small constant

hosts = as.data.frame(log10(metaphlan[, sele] + sc))    # log transformation to get normally distributed data


result.glm3 = foreach(i = 1:ncol(crass3), .combine = rbind) %:%

  foreach(j = 1:ncol(pheno2), .combine = rbind) %do% {

    if (colnames(pheno2)[j] == "antrop_age") {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, "antrop_age"] + pheno2[, "antrop_gender.F1M2"] + ., data = hosts, family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    } else if (colnames(pheno2)[j] == "antrop_gender.F1M2") {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, "antrop_age"] + pheno2[, "antrop_gender.F1M2"] + ., data = hosts, family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[3,1],
            SE = c1$coef[3,2],
            Z = c1$coef[3,3],
            P = c1$coef[3,4]
        )

    } else {

        c1 = summary(glm(crass3[, i] > 0 ~ pheno2[, j] + pheno2[, "antrop_age"] + pheno2[, "antrop_gender.F1M2"] + ., data = hosts, family = "binomial"))

        data.frame(
            pheno = colnames(pheno2)[j],
            taxon = colnames(crass3)[i],
            Beta = c1$coef[2,1],
            SE = c1$coef[2,2],
            Z = c1$coef[2,3],
            P = c1$coef[2,4]
        )

    }
    
  }

warnings()




# -------------------- multiple testing correction --------------------
result.quant$FDR = p.adjust(result.quant$P, method = "BH")
result.bin$FDR = p.adjust(result.bin$P, method = "BH")
result.glm1$FDR = p.adjust(result.glm1$P, method = "BH")
result.glm2$FDR = p.adjust(result.glm2$P, method = "BH")
result.glm3$FDR = p.adjust(result.glm3$P, method = "BH")




# -------------------- combine results --------------------
combined = data.frame(
    pheno = result.quant$pheno,
    taxon = result.quant$taxon,

    Rsp.quant = result.quant$Rsp,
    P.quant   = result.quant$P,
    FDR.quant = result.quant$FDR,

    Rsp.bin = result.bin$Rsp,
    P.bin   = result.bin$P,
    FDR.bin = result.bin$FDR,

    Beta.glm1 = result.glm1$Beta,
    SE.glm1   = result.glm1$SE,
    Z.glm1    = result.glm1$Z,
    P.glm1    = result.glm1$P,
    FDR.glm1  = result.glm1$FDR,

    Beta.glm2 = result.glm2$Beta,
    SE.glm2   = result.glm2$SE,
    Z.glm2    = result.glm2$Z,
    P.glm2    = result.glm2$P,
    FDR.glm2  = result.glm2$FDR,

    Beta.glm3 = result.glm3$Beta,
    SE.glm3   = result.glm3$SE,
    Z.glm3    = result.glm3$Z,
    P.glm3    = result.glm3$P,
    FDR.glm3  = result.glm3$FDR
)

combined = combined[order(combined$FDR.bin), ]




# -------------------- write table --------------------
combined$taxon <- sub("^o__", "", combined$taxon)
combined$taxon <- gsub("\\|[a-z]__", "|", combined$taxon)

write.table(combined, file = "LLD_pheno_assoc.txt", sep = "\t", row.names = F, quote = F)
