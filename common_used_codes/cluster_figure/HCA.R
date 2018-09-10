############make a volcano plot################################
###you should prepare a expression table
###the rownames is gene name and the colname is the samples name in the expression table
###you also should prepare a rt_sam_col table which contain three columns ("samples", "group", "color")
###UseMeathod: MakHCA(rt_hca, rt_sam_col = rt_sam_col, pnam, type = c("tiff", "pdf"), method = 'ward.D', cex = cex, width = NULL, height = NULL)
###20171018
###JunShang
###E-mail: shangjunv@163.com

##############test data#########################################
#rt_test <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch2/exp/MAQC_SampleAB_FPKM.txt",
#                      header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#rt_d <- rt_test[-which(apply(rt_test, 1, mean) == 0), ]
#rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
#colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
#
#rt_sam_col <- data.frame(cbind(colnames(rt_l), c("A", "A", "A", "B", "B", "B"), c("red", "red", "red", "blue", "blue", "blue")), 
#                         stringsAsFactors = FALSE)
#colnames(rt_sam_col) <- c("samples", "group", "color")
#MakHCA(rt_l, rt_sam_col = rt_sam_col, "MAQC", type = "tiff", method = 'ward.D', cex = 1.5, width = 1024, height = 716)
################################################################
library(base)
library(dendextend)
require(colorspace)

########################make a HCA plot#########################
MakHCA <- function(rt_hca, rt_sam_col = rt_sam_col, pnam, type = c("tiff", "pdf"), method = 'ward.D', cex = cex, width = NULL, height = NULL){
  #rt_hca is a expression table
  #the rownames is gene name and the colname is the samples name in the expression table
  #rt_sam_col table contain three columns ("samples", "group", "color")
  #pnam is used to choice your figure name and figure title
  
  rt_hclust <- hclust(d=dist(t(rt_hca)), method = method)
  rt_hclust_dend <- as.dendrogram(rt_hclust)
  new_order <- order.dendrogram(rt_hclust_dend)
  labels(rt_hclust_dend) <- colnames(rt_hca)[new_order]
  labels_colors(rt_hclust_dend) <- rt_sam_col$color[new_order]
  
  if(type == "tiff"){
    tiff(paste(pnam, "_HCA.tiff", sep = ""), width = width, height = height)
    par(cex = cex)
    plot(rt_hclust_dend, ylab = "Euclidean distance", main = paste("Clustered Samples", sep = ""))
    legend("topright", legend = unique(rt_sam_col$group), fill = unique(rt_sam_col$color), 
           text.col = unique(rt_sam_col$color))
    dev.off() 
  } else if (type == "pdf") {
    par(cex = cex)
    pdf(paste(pnam, "_HCA.pdf", sep = ""), width, height)
    plot(rt_hclust_dend, ylab = "Euclidean distance", main = paste("Clustered Samples", sep = ""))
    legend("topright", legend = unique(rt_sam_col$group), fill = unique(rt_sam_col$color), 
           text.col = unique(rt_sam_col$color))
    dev.off() 
  }
}
  
################################################################
