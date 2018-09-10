### when we have a gene list, we can merge all cancer expression data that contain the expression of these genes
### input file genes list and all cancers expression data


library(ggplot2)
library(plyr)
library(gplots)
library(plotly)
library(gmodels)  # the packages that PCA and boxplot need
library(gplots)  # heatmap
library(limma)  #normalizeQuantiles

gene_list <- read.table(file = "/Users/stead/Desktop/GDC/LUAD/ROC/tissue_special_gene/27_T_log2_special/gene_list.tsv", header = F, stringsAsFactors = F, sep = "\t")
gnam <- c(gene_list[, 1])
cnam <- c("LUAD", "BRCA", "BLCA", "ACC", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM",
                  "STAD", "TGCT", "THCA", "THYM", "UCEC", "LAML", "LGG")
k <- "0.."

splite_data <- function(rt, k){
  ### splite tumor and normal data and you can choose the output data base on k
  
  keepRow <- c(1: nrow(rt))            
  keepCol <- seq(1, ncol(rt)) 
  data <- rt[keepRow, keepCol]
  data <- log2(data + 0.001)
  dimnames <- list(rownames(data), colnames(data))
  data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames) 
  
  G.R <- sapply(strsplit(colnames(data), "\\."), "[", 4 ) # "["(a, 4) almost a[4]
  
  if(is.null(k)||k == ""){ 
    ### TCGA-xx-xxxx-xx-xxxx. The sample id is the condition choosing the data_out
    return(NULL)
  } else {
    data.index <- grep(k, G.R)
    data_out <- data[, data.index]
    return(data_out) 
  }
}

sort_data <- function(cnam, gnam, k){
  ### add cancer name as the group type
  rt <- read.table(file = paste(file = "/Users/stead/Desktop/GDC/", cnam, "/raw_data/", cnam, "_mRNA_FPKM.txt", sep = ""), header = T, na.strings = NA, row.names = 1, sep = "\t")
  rt_Tumor <- splite_data(rt, k)
  rt_Tumor_T <- rt_Tumor[gnam, ]
  data_T <- t(rt_Tumor_T)
  data_T_C <- cbind(c(rep(cnam, length(row.names(data_T)))), data_T)
  colnames(data_T_C)[1] <- "group"
  return(data_T_C)
}

data_T_A <- lapply(cnam, sort_data, gnam, k)
data_mat <- do.call("rbind", data_T_A)
data_mat_A <- cbind(data_mat[, 1], data.frame(apply(data_mat[, -1], 2, as.numeric)), stringsAsFactors = F)
colnames(data_mat_A)[1] <- "group"
write.table(data_mat_A, file = "cancer_T_genes_merge.txt", col.names = T, row.names = T, sep = "\t")


## dimnames <- list(rownames(data), colnames(data))
data_mat_A <- read.table(file = "/Users/stead/Desktop/GDC/LUAD/merge_result/27_log2exp_pac_heatmap/cancer_T_genes_merge.txt", header = T, na.strings = NA, row.names = 1, sep = "\t")
group <- data_mat_A[, 1]

### we must define the colors when you have many group
color_group <- c("mediumpurple4", "green", "yellow", "cyan3", "mediumpurple4", "yellow", 
                 "yellow", "cyan3", "cyan3", "cyan3", "mediumpurple4", "yellow", 
                 "green", "red", "blue", "cyan3", "mediumpurple4", "yellow", 
                 "cyan3", "yellow", "green", "cyan3", "mediumpurple4", "yellow", 
                 "cyan3", "mediumpurple4", "green")

### draw PCA 
### data_mat_A, the colums is the sample id, and the rows is the genes id 
data.a <- as.matrix(data_mat_A[, -1])
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame

### group can cut you can choose define. we can use them difine colors and shape 
pc$group <- group
#pc$cut <- cut 
xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file="pca_nq.tiff", width = 3600, height = 1800)
ggplot(pc, aes(PC1, PC2), main) +
  geom_point(size=15, aes(color = group)) + ###shape = cut, 
  scale_color_manual(values= color_group) +
  labs(x = xlab, y = ylab, title = "LUAD_SEG PCA") +
  theme_bw()+
  theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 50, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 50, hjust = 2, vjust = 3, face = 'bold'),
        legend.position = c(1, 1), legend.justification=c(1, 1), legend.background = element_rect(fill = 'white', colour = 'black', size = 1), 
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 3),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 3)) 

dev.off()

#data.G[1,][c(data.G[1,] == "Tumor")] <- 'darkblue'
#data.G[1,][c(data.G[1,] == "Normal")] <- 'green'
heatmap.cols <- c("red", "grey50", "grey51", "grey52", "grey53","grey54", "grey55", 
                  "grey56", "grey57", "grey58", "grey59", "grey60", "grey61", "grey62", "grey63",
                  "grey64", "grey65", "grey66", "grey67", "grey68",
                  "grey69", "grey70", "grey71", "grey72", "grey73", "grey74", "grey75")

cor_cut <- factor(group, labels = heatmap.cols)
cor_cut <- as.character(cor_cut)
data_mat_s <-t(data_mat_A[ ,-1])
data_top <- data_mat_s
tiff(file = "heatmap_DEG.tiff",width = 3600, height = 1800)
par(cex.lab = 3, cex.axis = 2, cex.main = 3 )
heatmap.2(data_top, col = colorpanel(99, "blue", "black", "red"), dendrogram = "column", 
          hclustfun = function(x){hclust(x, method = 'ward.D2')}, key=T, symkey=F, key.xlab = NA,
          key.title = "Normalied expres sion", keysize = 1, Colv = F, Rowv = T, trace = "none", density.info = "none",
          labRow = NA, labCol = NA, cexCol = 0.5, ColSideColors = cor_cut, scale = "row", na.rm = T,
          margin = c(14,10), colCol = heatmap.cols, srtCol = 90) 
legend("left", legend = cnam, pch = 15, col = heatmap.cols, cex = 3, text.col = "black") 
dev.off() 

### draw box plot
draw_box <- function(k, data_out){
  filename <- colnames(data_out)[k]
  tiff(file = paste(filename, ".tiff", sep = ""), width = 3600, height = 1800)
  p = ggplot(data = data_out, aes(group, data_out[,k])) +
    xlab("cancer_type") +  ylab("gene expression (log2FPKM)") +
    geom_boxplot(color="gray30", fill= color_group) +
    theme_bw()+
    theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), 
          axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          axis.text.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90),
          axis.title.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
          axis.text.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
          panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5))
  print(p)
  dev.off()   
} 

k <- c(2 :length(colnames(data_mat_A)))
lapply(k, draw_box, data_mat_A)

