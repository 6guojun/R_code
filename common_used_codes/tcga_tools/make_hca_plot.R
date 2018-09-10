library(dendextend)
library(dynamicTreeCut)
library(colorspace)

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

rt_hclust <- hclust(dist(t(rt_m_stageIA)), method = "ward.D2")
rt_hclust_dend <- as.dendrogram(rt_hclust)
new_order <- order.dendrogram(rt_hclust_dend)
rt_sam_col_d <- data.frame(cbind(colnames(rt_m_stageIA), c(rep("blue", length(grep("Normal", colnames(rt_m_stageIA)))), 
                                                rep("red", length(grep("Tumor", colnames(rt_m_stageIA)))))), stringsAsFactors = FALSE)

cols <- rt_sam_col_d$X2

pdf("hcluster.pdf", 20, 10)
plot(rt_hclust_dend)
colored_bars(cols, y_scale = 10, y_shift = 10)
dev.off()
