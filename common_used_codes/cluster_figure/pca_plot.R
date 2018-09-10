############make a volcano plot################################
###gplots package required
###you should prepare a expression table
###the rownames is gene name and the colname is the samples name in the expression table
###UseMeathod: MakPCA(rt_pca, rt_sam = rt_sam, nam = nam, ThemePCA, type = c("tiff", "pdf"), width = NULL, height = NULL)
###20170907
###JunShang
###E-mail: shangjunv@163.com

##############test data#########################################
#source("/Users/stead/Documents/SourceTree/R/RNA_seq/DEG_analysis/DEG_Test.R")
#rt <- read.table(file ="WJC_H9_FPKM.txt")
#rt <- rt[apply(rt, 1, mean) > 0, ]
#rt_sam <- read.table(file = "group.txt", header = TRUE, sep = "\t", row.name = NULL, stringsAsFactors = FALSE)
#DEGExp <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "DEGExp")
################################################################

library(gmodels)
########################make a PCA plot#########################
MakPCA <- function(rt_pca, rt_sam = rt_sam, nam = nam, ThemePCA, type = c("tiff", "pdf"), width = NULL, height = NULL){
  # rt_pca is the expresssion table 
  # rt_sam is the samples and group table
  # nam is used to set your figure name
  # you can choice which type of fiugre based on the parameter of type

  data.a <- t(as.matrix(rt_pca))  
  data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
  a <- summary(data.pca) 
  tmp <- a$importance  # a include 4 sections which contain importance 
  pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
  pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
  pc <- as.data.frame(a$x)  # convert to data.frame
  num_pos <- match(row.names(pc), rt_sam$samples, nomatch = NA)
  pc$group <-  rt_sam$group[num_pos]
  pc$color <-  rt_sam$color[num_pos]
  pc$group <- factor(pc$group, levels = unique(pc$group))
  print(c(pc$group, pc$color)) 
 
  xlab <- paste("PC1(", pro1, "%)", sep = "")
  ylab <- paste("PC2(", pro2, "%)", sep = "")
  if(type == "tiff"){
    tiff(file = paste(nam, "_PCA.tiff", sep = ""), width = width, height = height)
    PCA <- ggplot(pc, aes(PC1, PC2, color = group)) + 
      geom_point(size = 5) + scale_colour_manual(values = unique(pc$color)) +
      labs(x = xlab, y = ylab, title = "PCA") +
      ThemePCA
    print(PCA)
    dev.off()
  } else if (type == "pdf") {
    pdf(file = paste(nam, "PCA.pdf", sep = ""), width, height)
    PCA <- ggplot(pc, aes(PC1, PC2, color = group)) + 
      geom_point(size = 5) + scale_colour_manual(values = unique(pc$color)) +
      labs(x = xlab, y = ylab, title = "PCA") +
      ThemePCA
    print(PCA)
    dev.off()
  }
}
################################################################
