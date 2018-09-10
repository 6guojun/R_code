############make a heatmap plot################################
###gplots package required
###you should prepare a expression table
###the rownames is gene name and the colname is the samples name in the expression table
###you shoud also prepare GP1 and GP2 to show your group
###UseMeathod: MakHeatmap(rt_heat, rt_sam = NULL, GP1, GP2, type = c("tiff", "pdf"), width = NULL, height = NULL)
###20170907
###JunShang
###E-mail: shangjunv@163.com

library(gplots)
##############test data#########################################
#source("/Users/stead/Documents/SourceTree/R/RNA_seq/DEG_analysis/DEG_Test.R")
#rt <- read.table(file ="WJC_H9_FPKM.txt")
#rt <- rt[apply(rt, 1, mean) > 0, ]
#rt_sam <- read.table(file = "group.txt", header = TRUE, sep = "\t", row.name = NULL, stringsAsFactors = FALSE)
#DEGExp <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "DEGExp")
################################################################

#####################Make a Heatmap plot#######################
MakHeatmap <- function(rt_heat, rt_sam = NULL, GP1, GP2, type = c("tiff", "pdf"), width = NULL, height = NULL){
  # rt_heat is the expresssion table 
  # rt_sam is the samples and group table 
  # GP1 and GP2 reprsent names of two group
  # colnames of rt_heat can not have the same name of GP1 or GP2
  # you can choice which type figure based on the parameter of type
  
  if(!is.null(rt_sam)){
    num_pos <- match(colnames(rt_heat), rt_sam$samples, nomatch = NA)
    colnames(rt_heat) <- paste(colnames(rt_heat), rt_sam$group[num_pos], sep = "_")
    print("rt_sam is not NULL")
  } else if (is.null(rt_sam)) {
    print("rt_sam is NULL")
  }
  rt_heat_A <- apply(rt_heat, 2, FUN = function(x){(x-median(x))/sd(x)})
  colnames(rt_heat_A)[grep(GP1, colnames(rt_heat_A))] <- 'red'
  colnames(rt_heat_A)[grep(GP2, colnames(rt_heat_A))] <- 'blue'
  heatmap.cols <- colnames(rt_heat_A)
  colnames(rt_heat_A) <- colnames(rt_heat)
  if(type == "tiff"){
    tiff(file = paste(GP1, "_", GP2, "_Heatmap", ".tiff", sep = ""), width = width, height = height)
    par(mar = c(30, 10, 1, 60), cex.lab = 3, cex.axis = 2, cex.main = 3 )
    heatmap.2(rt_heat_A, col = colorpanel(99, "blue", "black", "red"), dendrogram = "both", 
              hclustfun = function(x){hclust(x, method = 'ward.D2')}, key = T, symkey = F, key.xlab = NA,
              key.title = "Normalied expression", keysize = 1, Colv = T, Rowv = T, trace = "none", density.info = "none",
              cexCol = 2, ColSideColors = heatmap.cols, scale = "row", labRow = NA, labCol = NA,
              margin = c(14,10), colCol = heatmap.cols, srtCol = 90) 
    legend("topright", legend = c(GP1, GP2), pch = 15, bty='n', col = c('red', 'blue'), cex = 2, text.col = "black") 
    dev.off() 
  } else if (type == "pdf") {
    pdf(file = paste(GP1, "_", GP2, "_Heatmap", ".pdf", sep = ""))
    heatmap.2(rt_heat_A, col = colorpanel(99, "blue", "black", "red"), dendrogram = "column", 
              hclustfun = function(x){hclust(x, method = 'ward.D2')}, key=T, symkey=F, key.xlab = NA,
              key.title = "Normalied expression", keysize = 1, Colv = T, Rowv = T, trace = "none", density.info = "none",
              labRow = NA, labCol = NA, cexCol = 2, ColSideColors = heatmap.cols, scale = "row",
              margin = c(14,10), colCol = heatmap.cols, srtCol = 90) 
    legend("topright", legend = c(GP1, GP2), pch = 15, bty='n', col = c('red', 'blue'), cex = 2, text.col = "black") 
    dev.off() 
  }
}
###############################################################
