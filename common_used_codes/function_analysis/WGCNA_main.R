library(WGCNA)

setwd(paste("/Users/stead/Desktop/subtype_analysis/signature/", cancer, sep = ''))
rt_exp <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt_T <- split_tcga_tn(rt_exp[, -1], sam_type = "tumor")



data_exp <- t(rt_GSEA_order_sym)
######select beta value######
powers1 <- c(seq(1, 10, by = 1))
RpowerTable <- pickSoftThreshold(data_exp, powerVector = powers1)[[2]]

pdf(file= "softThresholding.pdf")
par(mfrow = c(1, 2))
plot(RpowerTable[, 1], -sign(RpowerTable[, 3])*RpowerTable[, 2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[, 1], -sign(RpowerTable[, 3])*RpowerTable[, 2], labels=powers1,cex = 0.7,col = "red")
abline(h = 0.9, col = "red")
# this line corresponds to using an R^2 cut-off of h
plot(RpowerTable[, 1], RpowerTable[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n")
text(RpowerTable[, 1], RpowerTable[, 5], labels = powers1, cex = 0.7, col = "red")
dev.off()

for(beta1 in 1:10){
  ADJ = adjacency(data_exp, power = beta1)
  vis = exportNetworkToCytoscape(ADJ, edgeFile =paste('edge_', beta1, '.txt', sep = ''), nodeFile = paste('node_', beta1, '.txt', sep = ''), threshold = 0.9, includeColNames = TRUE)
  rt_node <- vis$edgeData
  write.table(rt_node, file = paste('related_edge_', beta1, '.txt', sep = '' ), col.names = TRUE, row.names = FALSE, sep = "\t")
  
}
