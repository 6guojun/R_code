############make a volcano plot################################
###gplots package required
###you should prepare a expression table
###the rownames is gene name and the colname is the samples name in the expression table
###UseMeathod: MakPCA(rt_pca, rt_sam = rt_sam, nam = nam, ThemePCA, type = c("tiff", "pdf"), width = NULL, height = NULL)
###20180911
###JunShang
###E-mail: shangjunv@163.com
library(gmodels)

ConSamMat <- function(rt_exp, grp_nam, len_a = len_a, len_b = len_b, col_a = col_a, col_b = col_b, g_cutoff = NULL){
  if(is.null(g_cutoff)){
    rt_exp$group[rt_exp[, grp_nam] > median(rt_exp[, grp_nam])] <- len_b
    rt_exp$group[rt_exp[, grp_nam] <= median(rt_exp[, grp_nam])] <- len_a
  } else {
    rt_exp$group[rt_exp[, grp_nam] > g_cutoff] <- len_b
    rt_exp$group[rt_exp[, grp_nam] <= g_cutoff] <- len_a
  }
  
  rt_exp$color <- NA
  rt_exp$color[grep(len_a, rt_exp$group)] <- col_a
  rt_exp$color[grep(len_b, rt_exp$group)] <- col_b
  
  rt_sam <- data.frame(cbind(row.names(rt_exp), rt_exp[, c('group', 'color')]), stringsAsFactors = FALSE)
  colnames(rt_sam) <- c('samples', 'group', 'color')
  
  return(rt_sam)
}


########################make a PCA plot#########################
MakPCA <- function(rt_pca, rt_sam = rt_sam, target_genes = NULL, nam = nam, phi = 40, pch = 18, len_a = len_a, len_b = len_b, 
                   len_position = 'bottom', width = NULL, height = NULL){
  # rt_pca is the expresssion table 
  # rt_sam is the samples and group table
  # nam is used to set your figure name
  # you can choice which type of fiugre based on the parameter of type
  if(is.null(target_genes)){
    rt_pca <- rt_pca
  }else if (!is.null(target_genes)){
    rt_pca <- rt_pca[, target_genes]
  }
  
  data.a <- as.matrix(rt_pca)
  data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
  a <- summary(data.pca) 
  tmp <- a$importance  # a include 4 sections which contain importance 
  pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
  pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100
  pro3 <- as.numeric(sprintf("%.3f", tmp[2,3]))*100# fetch the proportion of PC1 and PC2
  pc <- as.data.frame(a$x)  # convert to data.frame
  num_pos <- match(row.names(pc), rt_sam$samples_id, nomatch = NA)

  pc$group <-  rt_sam$group[num_pos]
  pc$color <-  rt_sam$color[num_pos]
  pc$group <- factor(pc$group, levels = unique(pc$group))
  
  xlab <- paste("PC1(", pro1, "%)", sep = "")
  ylab <- paste("PC2(", pro2, "%)", sep = "")
  zlab <- paste("PC3(", pro3, "%)", sep = "")
  
  pdf(file = paste(nam, "PCA.pdf", sep = ""), width, height)
  PCA <- scatter3D(pc$PC1, pc$PC2, pc$PC3,  xlab = xlab, ylab = ylab, zlab = zlab, bty = 'g', cex = 1.5,
                   pch  = pch, phi = phi, col = as.character(unique(pc$color)),  colvar = as.integer(pc$group), colkey = FALSE) 
    legend(len_position, legend = levels(pc$group), col = as.character(unique(pc$color)), pch = 18, inset = -0.17, xpd = TRUE, horiz = TRUE, cex = 1.5)
  print(PCA)
  dev.off()
}
################################################################
