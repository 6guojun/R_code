### input file miRNA "Ttest_DEG_EXP", miRNA_id_transform.txt, mRNA Ttest_DEG_EXP.txt
### 

miRNA_DEG_exp <- read.table(file = "/Users/stead/Desktop/GDC/LUSC/DEG_miRNA/Ttest_DEG_EXP.txt", header = T, row.names = NULL, sep = "\t")
miRNA_id_T <- read.table(file = "/Users/stead/Desktop/miRNA_id_transform.txt", header = T, sep = "\t")

splite_mi_data <- function(miRNA_DEG, miRNA_id, t){
  if(t == "Tumor"){
    index <- grep('Tumor', colnames(miRNA_DEG_exp)) ## get the tumor group index
    miRNA_DEG_exp_T <- miRNA_DEG_exp[, index]
    miRNA_DEG_exp_T <- cbind(miRNA_DEG_exp[, 1], miRNA_DEG_exp_T)
    
    miRNA_DEG_id_T <- merge(miRNA_id_T, miRNA_DEG_exp_T, by = 1)
    miRNA_DEG_nam <- miRNA_DEG_id_T[, c(-1, -(3:6))]
    row.names(miRNA_DEG_nam) <- miRNA_DEG_nam[, 1]
    miRNA_DEG_T_out <- data.frame(apply(miRNA_DEG_nam[, -1], 2, as.numeric), stringsAsFactors = F)
    row.names(miRNA_DEG_T_out) <- miRNA_DEG_nam[, 1]
    return(miRNA_DEG_T_out)
  } else{
    index <- grep('Normal',colnames(miRNA_DEG_exp))
    miRNA_DEG_exp_T <- miRNA_DEG_exp[, index]
    miRNA_DEG_exp_T <- cbind(miRNA_DEG_exp[, 1], miRNA_DEG_exp_T)
    
    miRNA_DEG_id_T <- merge(miRNA_id_T, miRNA_DEG_exp_T, by = 1)
    miRNA_DEG_nam <- miRNA_DEG_id_T[, c(-1, -(3:6))]
    row.names(miRNA_DEG_nam) <- miRNA_DEG_nam[, 1]
    miRNA_DEG_T_out <- data.frame(apply(miRNA_DEG_nam[, -1], 2, as.numeric), stringsAsFactors = F)
    row.names(miRNA_DEG_T_out) <- miRNA_DEG_nam[, 1]
    return(miRNA_DEG_T_out)
  }
}
miRNA_DEG_T_out <- splite_mi_data(miRNA_DEG_exp, miRNA_id_T, "Tumor")



mRNA_DEG_exp <- read.table(file = "/Users/stead/Desktop/GDC/LUSC/DEG/Ttest_DEG_EXP.txt", header = T, row.names = NULL, sep = "\t", stringsAsFactors = F) 
#rt_T <- read.table(file = "/Users/stead/Desktop/GDC/LUSC/ROC/tissue_special_gene/27_T_log2_special/gene_list.tsv", header = F, row.names = NULL, sep = "\t", stringsAsFactors = F)
rt_T <- read.table(file = "/Users/stead/Desktop/GDC/LUSC/DEG/DEGs_id_symbol_exp.txt", header = T, row.names = 1,  sep = "\t", stringsAsFactors = F)
rt_T_nam <- data.frame(cbind(rt_T$V1, rt_T$hgnc_symbol), stringsAsFactors = F)

mRNA_DEG_T_exp <- merge(rt_T_nam , mRNA_DEG_exp, by = 1)
mRNA_DEG_T_exp <- mRNA_DEG_T_exp[-c(which(mRNA_DEG_T_exp[, 2] == "")), ]
mRNA_DEG_T_T <- data.frame(apply(mRNA_DEG_T_exp[, c(-1 :-2)], 2, as.numeric))
row.names(mRNA_DEG_T_T) <- mRNA_DEG_T_exp[, 2]


tumor.index <- grep('Tumor', colnames(mRNA_DEG_T_T)) ## get the tumor group index
tumor.index
normal.index <- grep('Normal',colnames(mRNA_DEG_T_T))
normal.index


mRNA_DEG_T_out <- mRNA_DEG_T_T[, tumor.index]

colnames(mRNA_DEG_T_out) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*", "\\1\\.\\2\\.\\3\\.\\4", colnames(mRNA_DEG_T_out))
colnames(mRNA_DEG_T_out) <- substring(colnames(mRNA_DEG_T_out), 1, 15)
colnames(miRNA_DEG_T_out) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*", "\\1\\.\\2\\.\\3\\.\\4", colnames(miRNA_DEG_T_out))
colnames(miRNA_DEG_T_out) <- substring(colnames(miRNA_DEG_T_out), 1, 15)

sameSample <- intersect(colnames(mRNA_DEG_T_out), colnames(miRNA_DEG_T_out))
merge <- rbind(id = sameSample, miRNA_DEG_T_out[, sameSample], mRNA_DEG_T_out[, sameSample])
merge <- merge[-1, ]
write.table(merge, file = "merge.txt", sep="\t", quote = F, col.names = F)

mRNALabel <- cbind(rownames(mRNA_DEG_T_out),"mRNA")
miRNALabel <- cbind(rownames(miRNA_DEG_T_out),"miRNA")
nodeLabel <- rbind(c("ID","Classify"), miRNALabel, mRNALabel)
write.table(nodeLabel, file="nodeLabel.txt", sep="\t", quote=F, col.names=F, row.names=F)


library(WGCNA)
library(doParallel)

data_merge <- data.frame(apply(merge, 2, as.numeric), stringsAsFactors = F)
row.names(data_merge) <- row.names(merge)
datExpr <- t(data_merge)

######select beta value######
powers1 <- c(seq(1, 10, by=1), seq(12, 30, by=2))
RpowerTable <- pickSoftThreshold(datExpr, powerVector=powers1)[[2]]
cex1 <- 0.6
pdf(file="softThresholding.pdf")
par(mfrow=c(1, 2))
plot(RpowerTable[,1], -sign(RpowerTable[, 3])*RpowerTable[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n")
text(RpowerTable[,1], -sign(RpowerTable[, 3])*RpowerTable[, 2], labels = powers1,cex=cex1, col = "red")
### this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red")
dev.off()

beta1 <- 5
ADJ <- adjacency(datExpr,power = beta1)
vis <- exportNetworkToCytoscape(ADJ, edgeFile = "edge.txt", nodeFile = "node.txt", threshold = 0.5)




