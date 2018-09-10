library(survival)
library(ggplot2)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(annotables)
library(pathview)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/survival_analysis")
cancer="LUAD"

rt <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                 sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt$Ensembl_ID <- matrix(unlist(strsplit(rt$Ensembl_ID, "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(rt$Ensembl_ID, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
rt_sym <- merge(eg, rt, by = 1)
row.names(rt_sym) <- paste(rt_sym$ENSEMBL, rt_sym$SYMBOL, sep = "_")
rt_sym <- rt_sym[, -c(1, 2)]

source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")
rt_T <- split_tcga_tn(rt_sym, sam_type = "tumor", split_type = "[.]")

sam_G <- matrix(unlist(strsplit(colnames(rt_sym), "[.]")), ncol = 4, byrow = TRUE)[, 4]
pos_T <- c(grep("01", sam_G), grep("02", sam_G))


##import clinical data
rt_sur <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/phenotype/TCGA-", cancer, ".survival.tsv", sep = ""),
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt_sur_old <- read.table(file = "/Users/stead/Desktop/PGX_group/lungcancer_project/20161124LUAD_RNA_seq/20161124LUAD_FPKM/survival/clinical_data_survival.csv",
                         sep = ",", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)

sam_sur <- matrix(unlist(strsplit(rt_sur_old$id, "[.]")), ncol = 4, byrow = TRUE)
rt_sur_old$id <- paste(sam_sur[, 1], sam_sur[, 2], sam_sur[, 3], sep = "-")
pari_sam <- match(colnames(rt_T), rt_sur_old$id , nomatch = 0)
pari_sam_r <- which(pari_sam != 0)
print(pari_sam_r)

rt_T_M <- rt_T[, pari_sam_r]
rt_sur_M <- rt_sur_old[pari_sam, ]

###draw survival figure
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/survival_analysis_2.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/survival_analysis.R")


GetSur <- function(gnam, rt_T_M){
  
  Rt_Exp_CLi <- data.frame(cbind(rt_sur_M[, c("status", "time")], t(rt_T_M[gnam, ])), stringsAsFactors = FALSE)
  colnames(Rt_Exp_CLi) <- c("OS_Status", "OS_Time", gnam)
  SurMainFun(gnam, Rt_Exp_CLi, DatType = "ConType", feature = "OS", OutType = "SurF", SurvType = "ALL")
  Pvalue <- SurMainFun(gnam, Rt_Exp_CLi, DatType = "ConType", feature = "OS", OutType = "SurP", SurvType = "ALL")
  return(Pvalue)
}

gnam <- row.names(rt_T_M)[1:100]
pval_list <- lapply(gnam, GetSur, rt_T_M)
rt_pval <- do.call(rbind, pval_list)

rt_pval_0.05 <- data.frame(rt_pval[is.na(rt_pval[, 2]), ], stringsAsFactors = FALSE)[, c(1, 3)]
rt_pval_0.05[, 2] <- as.numeric(rt_pval_0.05[, 2])
colnames(rt_pval_0.05) <- c("gene_id", "pvalue")
rt_pval_0.01 <- rt_pval_0.05[rt_pval_0.05$pvalue < 0.01, ]
rt_pval_0.001 <- rt_pval_0.05[rt_pval_0.05$pvalue < 0.001, ]
write.table(rt_pval_0.05, file = paste(cancer, "_P_005.txt", sep = ""), col.names = TRUE, row.names = FALSE)
write.table(rt_pval_0.01, file = paste(cancer, "_P_001.txt", sep = ""), col.names = TRUE, row.names = FALSE)
write.table(rt_pval_0.001, file = paste(cancer, "_P_0001.txt", sep = ""), col.names = TRUE, row.names = FALSE)
