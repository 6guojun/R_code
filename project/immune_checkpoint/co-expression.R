library(WGCNA)
setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/co_expression")
rt <- read.table("/Users/stead/Desktop/PD-L1_and_TMI_type/UVM/RNA_seq/TCGA-UVM.htseq_fpkm.tsv",sep="\t", row.names=1, header=T, stringsAsFactors = FALSE)
row.names(rt)<- matrix(unlist(strsplit(row.names(rt), "[.]")), ncol = 2, byrow = TRUE)[, 1]
rt_d <- rt[apply(rt, 1, function(x){median(x) > 0}), ]
datExpr <- t(rt_d[1:1000, ])

######select beta value######
powers1 <- c(seq(1, 10, by = 1))
RpowerTable <- pickSoftThreshold(datExpr, powerVector=powers1)[[2]]

pdf(file= paste(cancer, "_softThresholding.pdf", sep = ""))
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers1,cex = 0.7,col = "red")
abline(h = 0.9, col = "red")
# this line corresponds to using an R^2 cut-off of h
plot(RpowerTable[, 1], RpowerTable[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n")
text(RpowerTable[, 1], RpowerTable[, 5], labels = powers1, cex = 0.7, col = "red")
dev.off()


beta1 = 6
ADJ = adjacency(datExpr, power = beta1)

##threshold 0.9
vis_0.9 = exportNetworkToCytoscape(ADJ,edgeFile = paste(cancer, "_edge_0.9.txt", sep = ""),
                             nodeFile =  paste(cancer, "_node_0.9.txt", sep = ""), threshold = 0.9, includeColNames = TRUE)
edge_0.9 <- vis_0.9$edgeData
edge_0.9_gene <- edge_0.9[c(grep("ENSG00000188389", edge_0.9$fromNode), grep("ENSG00000188389", edge_0.9$toNode)), ]
write.table(edge_0.9_gene, file = paste(cancer,  "_edge_0.9_PD1.txt", sep = ""), col.names = TRUE, row.names = FALSE)


##threshold 0.8
vis_0.8 = exportNetworkToCytoscape(ADJ,edgeFile = paste(cancer,  "_edge_0.8.txt", sep = ""),
                               nodeFile =  paste(cancer,  "_node_0.8.txt", sep = ""), threshold = 0.8, includeColNames = TRUE)
edge_0.8 <- vis_0.8$edgeData
edge_0.8_gene <- edge_0.8[c(grep("ENSG00000188389", edge_0.8$fromNode), grep("ENSG00000188389", edge_0.8$toNode)), ]
write.table(edge_0.8_gene, file = paste(cancer,  "_edge_0.8_PD1.txt", sep = ""), col.names = TRUE, row.names = FALSE)

##threshold 0.7
vis_0.7 = exportNetworkToCytoscape(ADJ,edgeFile = paste(cancer, "_edge_0.7.txt", sep = ""),
                               nodeFile =  paste(cancer, "_node_0.7.txt", sep = ""), threshold = 0.7, includeColNames = TRUE)
edge_0.7 <- vis_0.7$edgeData
edge_0.7_gene <- edge_0.7[c(grep("ENSG00000188389", edge_0.7$fromNode), grep("ENSG00000188389", edge_0.7$toNode)), ]
write.table(edge_0.7_gene, file = paste(cancer,  "_edge_0.7_PD1.txt", sep = ""), col.names = TRUE, row.names = FALSE)

##threshold 0.6
vis_0.6 = exportNetworkToCytoscape(ADJ,edgeFile = paste(cancer, "_edge_0.6.txt", sep = ""),
                               nodeFile =  paste(cancer, "_node_0.6.txt", sep = ""), threshold = 0.6, includeColNames = TRUE)
edge_0.6 <- vis_0.6$edgeData
edge_0.6_gene <- edge_0.6[c(grep("ENSG00000188389", edge_0.6$fromNode), grep("ENSG00000188389", edge_0.6$toNode)), ]
write.table(edge_0.6_gene, file = paste(cancer,  "_edge_0.6_PD1.txt", sep = ""), col.names = TRUE, row.names = FALSE)



