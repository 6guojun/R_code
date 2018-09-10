library(clusterProfiler)
library(ggplot2)
library(ggcorrplot)
library(ggdendro)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint")

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

count_cancer_cor <- function(cnam, ctr_gene, test_gene, data_cor, out_type = c("name", "noname")){
  cancer_pos <- grep(cnam, data_cor$cancer_type)
  ctr_gene_exp <- as.numeric(data_cor[, ctr_gene][cancer_pos])
  test_gene_exp <- as.numeric(data_cor[, test_gene][cancer_pos])
  cor_value <- round(cor(ctr_gene_exp, test_gene_exp), digits = 3)
  if(out_type == "noname"){
    return(cor_value)
  } else if (out_type == "name"){
    return(c(paste(cnam, ctr_gene, test_gene, sep = "_"), cor_value))
  }
}

#count_cor("GBM", "CTLA4", "CD80", rt_cor)

count_gene_cor <- function(gnam, cancer_list, ctr_gene, data_cor, out_type){
  cnam = cancer_list
  cor_value <- lapply(cnam, count_cancer_cor, ctr_gene, gnam, data_cor, out_type = out_type)
  cancer_gene_cor <- do.call(cbind, cor_value)
  return(cancer_gene_cor)
}


GeneToGenesCor <- function(ctrl_gene, cancer_list, rt_cor, out_type){
  gnams <- colnames(rt_cor)[-match(c("cancer_type", ctrl_gene), colnames(rt_cor))]
  ic_cor_list_nam <- lapply(gnams, count_gene_cor, cancer_list, ctrl_gene, rt_cor, "name")
  ic_cor_name <- do.call(rbind, ic_cor_list_nam)
  
  ic_cor_list_noname <- lapply(gnams, count_gene_cor, cancer_list, ctrl_gene, rt_cor, "noname")
  ic_cor_noname <- do.call(rbind, ic_cor_list_noname)
  colnames(ic_cor_noname) <- cancer_list
  row.names(ic_cor_noname)  <- gnams
  if(out_type == "noname"){
    return(ic_cor_noname)
  } else {
    return(ic_cor_name)
  }
}

theme_temp <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 10, colour = 'black', vjust = 0.5, hjust = 0.5, angle = 45), 
  axis.title.y = element_blank(),
  axis.text.y = element_text(size = 10, colour = 'black', vjust = 0.5, hjust = 0.5, angle = 1), 
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.text = element_text(colour = 'black', angle = 1, size = 10, hjust = 2, vjust = 3), 
  legend.position = "right", legend.direction = "vertical")


###do CTLA4 and immune point correaltion
ic_cor_noname_CTLA4 <- GeneToGenesCor("CTLA4", cancer_list, rt_cor, "noname")
#do cluster 
ord_CTLA4_dendro <- as.dendrogram(hclust(d = dist(x = t(ic_cor_noname_CTLA4))))
dendro_plot_CTLA4 <- ggdendrogram(data = ord_CTLA4_dendro, rotate = FALSE)
otter_order_CTLA4 <- order.dendrogram(ord_CTLA4_dendro)

#Re-order heatmap columns to match dendrogram
ic_cor_noname_CTLA4_ord <- ic_cor_noname_CTLA4[, otter_order_CTLA4]

#Putting clusetr and heatmap all together
melted_cormat_CTLA4 <- melt(ic_cor_noname_CTLA4_ord)
pdf(file = paste("CTLA4_icheckpoint_cluster", ".pdf", sep = ""), 10, 15)
ggheatmap_CTLA4 <- ggplot(melted_cormat_CTLA4, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "correlation") +
  geom_text(aes(Var2, Var1, label = value), color = "black",  size = 2) +
  theme_temp #coord_fixed ###fix cell size

grid.newpage()
print(ggheatmap_CTLA4, vp=viewport(0.8, 0.8, x = 0.45, y = 0.4)) #x change x positiion, y change y position, 
print(dendro_plot_CTLA4, vp=viewport(0.43, 0.2, x = 0.43, y = 0.9)) # 0.5change cluster plot wide
dev.off()


###do CD274 and immune point correaltion
ic_cor_noname_CD274 <- GeneToGenesCor("CD274", cancer_list, rt_cor, "noname")

#do cluster 
ord_CD274_dendro <- as.dendrogram(hclust(d = dist(x = t(ic_cor_noname_CD274))))
dendro_plot_CD274 <- ggdendrogram(data = ord_CD274_dendro, rotate = FALSE)
otter_order_CD274 <- order.dendrogram(ord_CD274_dendro)

#Re-order heatmap columns to match dendrogram
ic_cor_noname_CD274_ord <- ic_cor_noname_CD274[, otter_order_CD274]

#Putting clusetr and heatmap all together
melted_cormat_CD274 <- melt(ic_cor_noname_CD274_ord)
pdf(file = paste("CD274_icheckpoint_cluster", ".pdf", sep = ""), 10, 15)
ggheatmap_CD274 <- ggplot(melted_cormat_CD274, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "correlation") +
  geom_text(aes(Var2, Var1, label = value), color = "black",  size = 2) +
  theme_temp #coord_fixed ###fix cell size

grid.newpage()
print(ggheatmap_CD274, vp=viewport(0.8, 0.8, x = 0.45, y = 0.4)) #x change x positiion, y change y position, 
print(dendro_plot_CD274, vp=viewport(0.43, 0.2, x = 0.43, y = 0.9)) # 0.5change cluster plot wide
dev.off()

####do PDCD1 and immune point correaltion

ic_cor_noname_PDCD1 <- GeneToGenesCor("PDCD1", cancer_list, rt_cor, "noname")

#do cluster 
ord_PDCD1_dendro <- as.dendrogram(hclust(d = dist(x = t(ic_cor_noname_PDCD1))))
dendro_plot_PDCD1 <- ggdendrogram(data = ord_PDCD1_dendro, rotate = FALSE)
otter_order_PDCD1 <- order.dendrogram(ord_PDCD1_dendro)

#Re-order heatmap columns to match dendrogram
ic_cor_noname_PDCD1_ord <- ic_cor_noname_PDCD1[, otter_order_PDCD1]

#Putting clusetr and heatmap all together
melted_cormat_PDCD1 <- melt(ic_cor_noname_PDCD1_ord)
pdf(file = paste("PDCD1_icheckpoint_cluster", ".pdf", sep = ""), 10, 15)
ggheatmap_PDCD1 <- ggplot(melted_cormat_PDCD1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "correlation") +
  geom_text(aes(Var2, Var1, label = value), color = "black",  size = 2) +
  theme_temp #coord_fixed ###fix cell size

grid.newpage()
print(ggheatmap_PDCD1, vp=viewport(0.8, 0.8, x = 0.45, y = 0.4)) #x change x positiion, y change y position, 
print(dendro_plot_PDCD1, vp=viewport(0.43, 0.2, x = 0.43, y = 0.9)) # 0.5change cluster plot wide
dev.off()

