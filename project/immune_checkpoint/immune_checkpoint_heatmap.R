library(clusterProfiler)
library(ggplot2)
library(ggcorrplot)
library(ggdendro)
library(ComplexHeatmap)
library(circlize)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint")
color_type_table <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/color.txt", sep = "\t", stringsAsFactors = FALSE)
color_type <- color_type_table$V2
cancer_list <- c("GBM", "UVM", "LUAD")
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
#colnames(rt_caner_all) <- c("cancer_type", ensemble_gene)
rt_ic <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/immune_checkpoint.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ens_ic <- eg$ENSEMBL[match(rt_ic$genes, eg$SYMBOL, nomatch = 0)]
ic_exp <- rt_caner_all[, -1][, match(ens_ic, ensemble_gene, nomatch = 0)]
colnames(ic_exp) <- rt_ic$genes
rt_cor <- data.frame(rt_caner_all$cancer, ic_exp, stringsAsFactors = FALSE)
colnames(rt_cor)[1] <- "cancer_type"

cancer_type <- data.frame(type =  rt_cor[, 1], stringsAsFactors = FALSE)
col_type <- factor(cancer_type$type, levels = unique(cancer_type$type),  labels = color_type[1 : length(unique(cancer_type$type))])
mat <- apply(rt_cor[, -1], 2, as.numeric)

pdf(file = paste("test", "_immune_checkpoint.pdf", sep = ""), 10, 15)
ha = rowAnnotation(df = cancer_type, gp = gpar(col_type), width = unit(1, "cm"))
ht_list <- Heatmap(mat, name = "immune_checkpoint", col = colorRamp2(c(0, 4, 8), c( "cyan4", "white", "red")),
                        cluster_columns = TRUE, show_row_dend = TRUE, rect_gp = gpar(col= "white"), show_column_names = TRUE, show_column_dend = TRUE, 
                        cluster_rows = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 15), row_names_max_width = unit(18, "cm"),
                        column_title = 'number of immune pathway genes')
ht_all <- ha + ht_list
draw(ht_all, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()
