library(survival)
library(ggplot2)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(annotables)
library(pathview)
library(RColorBrewer)
library(scales)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/gene_all_cancer_exp")

cancer="GBM"
cancer_list <- c("GBM", "LUAD")
###read all tcga cancer expressino table and merge
###cancer list





source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")

GetBigMat <- function(cancer_list){
  cancer <- cancer_list[1]
  rt_on <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""), 
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  rt_T <- split_tcga_tn(rt_on, sam_type = "tumor", split_type = "[.]")
  rt_on_m <- data.frame(cbind(cancer, t(rt_T)), stringsAsFactors = FALSE)

  for(i in 2: length(cancer_list)){
    ###add the other cancer types
    cancer_add <- cancer_list[i]
    rt_add <-read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer_add, "/RNA_seq/TCGA-", cancer_add, ".htseq_fpkm.tsv", sep = ""), 
                        sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rt_add_T <- split_tcga_tn(rt_add, sam_type = "tumor", split_type = "[.]")
    rt_add_m <- data.frame(cbind(cancer_add, t(rt_add_T)), stringsAsFactors = FALSE)
    colnames(rt_add_m)[1] <- "cancer"
    if(all(colnames(rt_on_m) == colnames(rt_add_m))){
      rt_on_m <- data.frame(rbind(rt_on_m, rt_add_m), stringsAsFactors = FALSE)
    } else {
      stop('the gene names and order of all cancer must be same')
    }
  }
  return(rt_on_m)
}

rt_caner_all <- GetBigMat(cancer_list)

ensemble_gene <- matrix(unlist(strsplit(colnames(rt_caner_all)[-1], "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(ensemble_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
symbol_gene <- eg$SYMBOL[match(ensemble_gene, eg$ENSEMBL, nomatch = 0)]
colnames(rt_caner_all) <- c("group", paste(ensemble_gene, symbol_gene, sep = "_"))

###make a boxlot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
for(nam in colnames(rt_caner_all)[-1]){
  rt_box <- rt_caner_all[, c("group", nam)]
  MakBoxPlotOder(nam, "group", rt_box, width = 6, height = 6, theme_E) 
}



