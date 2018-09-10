###make a assess to MCPcount score and immune checkpoint correlationship
###
###shangjun 
###20180412
###E-mail: shangjunv@163.com
library(devtools)
library(MCPcounter)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggcorrplot)
library(ggdendro)
library(reshape2)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/MCP/each_cancer")
cancer="GBM"
gnam="PDCD1"
color_type_table <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/color.txt", sep = "\t", stringsAsFactors = FALSE)
color_type <- color_type_table$V2

####################################################
#get cancer exprssion 
#get the number of mutation and neo-ag
####################################################

####cancer expression and convert
MCPScoreMatrix <- function(cancer, gnam){
  rt <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                   sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  row.names(rt) <- matrix(unlist(strsplit(row.names(rt), "[.]")), ncol = 2, byrow = TRUE)[, 1]
  eg = bitr(row.names(rt), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  rt_sym <- merge(eg, cbind(row.names(rt), rt), by = 1)[, -1] 
  
  source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")
  rt_T <- split_tcga_tn(rt_sym[, -1], sam_type = "tumor", split_type = "[.]")
  #rt_N <- split_tcga_tn(rt_sym[, -1], sam_type = "normal")
  
  
  ###################################################
  #get score matrix
  #
  ###################################################
  make_mcp <- function(rt_mcp){
    ###your gene ID should be set gene sumbol
    load("/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/MCP_probesets.Rdata")  
    load("/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/MCP_genes.Rdata") 
    rt_score <- MCPcounter.estimate(rt_mcp, featuresType = c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[2],
                                    probesets = probesets, genes = genes)
    return(rt_score)
  }
  
  
  ###################################################
  #MCP counter all samples included tumor 
  #
  ###################################################
  mcp_gene <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/MCP/MCP_counter_signature.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  mcp_gene_pos <- match(mcp_gene$HUGO.symbols, rt_sym$SYMBOL, nomatch = 0)
  
  ###the rownames of rt_T/N is the same with rt_sym
  rt_T_cmp <- rt_T[mcp_gene_pos, ]
  row.names(rt_T_cmp) <- rt_sym$SYMBOL[mcp_gene_pos]
  
  rt_score_T <- make_mcp(rt_T_cmp)
  #tumor
  pdf(file = paste(cancer, "cmp_sore_tumor_raw.pdf", sep = "_"), 10, 5)
  p <- heatmap(as.matrix(rt_score_T),col=colorRampPalette(c("blue","white","red"))(100))
  print(p)
  dev.off()
  
  gexp  <- as.numeric(rt_T[which(rt_sym$SYMBOL == gnam), ])
  names(gexp) <- colnames(rt_T)
  rt_score_list <- list(rt_score_T, gexp)
  return(rt_score_list)
}


MCPGeneCor <- function(cancer, gnam){
  rt_score_list <- MCPScoreMatrix(cancer, gnam)
  rt_score_T <- rt_score_list[[1]]
  gexp <- rt_score_list[[2]]
  cor_value <- round(cor(gexp, t(rt_score_T)), digits = 2) 
  row.names(cor_value) <- cancer
  print(cancer)
  #make a correlation plot
  source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/Mut_correlation.R")
  rt_score_genes <- rbind(rt_score_T, gexp) 
  MakCorPlot(t(rt_score_genes), Tnam = paste(cancer, "_gene_score_cor", sep = ""), height = 15, width = 15, cex.axis = 2)
  return(cor_value)
}

###get MCP and gene expression correlation matrix
cancer <- c("UVM", "GBM", "OV")
cor_value_list <- lapply(cancer, MCPGeneCor, "PDCD1")
rt_cor_all <- do.call(rbind, cor_value_list)

rt_cor_all <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/PD1_all_cancer_MCP.txt", sep = "\t")
###get cluster matrix
source("/Users/stead/Documents/SourceTree/R/common_used_codes/get_cluster_mat.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/theme_cor_heatmap.R")
rt_cor_all_t <- t(rt_cor_all)
rt_reorder_list <- MakCluterMat(rt_cor_all_t, method = "complete", cluster_dir = "all")
rt_reorder <- rt_reorder_list$reorder_mat
ord_dendro_V <- rt_reorder_list$ord_dendro_V
ord_dendro_H <- rt_reorder_list$ord_dendro_H
dendro_plot_V <- ggdendrogram(data = ord_dendro_V, rotate = FALSE, labels = FALSE)
dendro_plot_H <- ggdendrogram(data = ord_dendro_H, rotate = TRUE, labels = FALSE)
melted_cormat <- melt(rt_reorder)

pdf(file = paste("icheckpoint_cluster", ".pdf", sep = ""), 40, 8)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "correlation") +
  geom_text(aes(Var2, Var1, label = value), color = "black",  size = 8) +
  theme_cor_heatmap #coord_fixed ###fix cell size
grid.newpage()
print(ggheatmap, vp=viewport(0.8, 0.8, x = 0.4, y = 0.4)) #x change x positiion, y change y position, 
print(dendro_plot_V, vp=viewport(w = 0.74, h = 0.15, x = 0.45, y = 0.9)) # 0.5change cluster plot wide
print(dendro_plot_H, vp=viewport(w = 0.05, h = 0.66, x = 0.83, y = 0.46))
dev.off()

