 ###sort data by FUSCC_RNAseq_exp_all.R

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")

library(dendextend)
library(dendextendRcpp)
library(dynamicTreeCut)
library(colorspace)

rt_l_clust <- rt_l_d

mgsub <- function(pattern, replacement, x, ...) {
  n = length(pattern)
  if (n != length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result = x
  for (i in 1:n) {
    result[grep(pattern[i], x, ...)] = replacement[i]
  }
  return(result)
}


rt_hclust <- hclust(dist(t(rt_l_clust)), method = "ward.D2")
rt_hclust_dend <- as.dendrogram(rt_hclust)
new_order <- order.dendrogram(rt_hclust_dend)

cols <- rt_sam_col_d$color
batch <- rt_sam_col_d$batch_s
color_set <- c("red", "cornflowerblue", "seagreen3", "purple3", "orange2", "yellow", "orange4", "violet", "grey", "wheat1", "gold",
             "cyan", "black", "white", "blue", "tomato")#, "plum2")
names(batch) <- mgsub(unique(batch), color_set, batch)


labels(rt_hclust_dend) <- colnames(rt_l_clust)[new_order]
labels_colors(rt_hclust_dend) <- rt_sam_col_d$color[new_order]

tiff(paste("28remove", "_HCA_bold.tiff", sep = ""), width = 2048, height = 1068)
par(mar = c(20, 4, 4, 2), cex = 1,  lwd = 2)
plot(rt_hclust_dend)
colored_bars(cbind(names(batch), cols), rt_hclust_dend, y_scale = 550, y_shift = -200)
legend("topright", legend = unique(rt_sam_col_d$group), fill = unique(rt_sam_col_d$color), 
       text.col = unique(rt_sam_col_d$color))
legend("right", legend = unique(rt_sam_col_d$batch_s), fill = unique(names(batch)), 
       text.col = unique(names(batch)))
dev.off()



