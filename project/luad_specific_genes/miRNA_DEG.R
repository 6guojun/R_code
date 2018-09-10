library(ggplot2)
library(plyr)
library(gplots)
library(plotly)
library(gmodels)  # the packages that PCA and boxplot need
library(gplots)  # heatmap
library(limma)  #normalizeQuantiles

result_dir = "/Users/stead/Desktop/GDC/LUSC/DEG_miRNA"
data_dir = "/Users/stead/Desktop/GDC/"

setwd(paste(data_dir, "LUSC/", "DEG_miRNA", sep = ""))
rt <- read.table(file = paste(data_dir, "LUSC", "/raw_data", "/LUSC_miRNA_HiSeq_gene", sep = ""), header = T, row.names = 1, na.strings = NA, stringsAsFactors = F, sep = "\t")
colnames(rt) <- colnames(rt)[-1] ###only first run
keepRow <- c(1: nrow(rt))            
keepCol <- seq(1, ncol(rt)) ## product the data of read per millions
data <- rt[keepRow, keepCol]
dimnames <- list(rownames(data), colnames(data))
data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames) 
#data <- ### convert to matrix
dim(data) 

G.R <- sapply(strsplit(colnames(data), "\\."), "[", 4 ) # "["(a, 4) almost a[4]
G.R <- sapply(strsplit(G.R, ""), "[", 1 ) 
G.R.R <- gsub("2", "1", G.R)    ## convert 2 to 1
G.R.R <- sub("0", "Tumor", G.R.R)  ## "T" replace "0"
G.R.R <- sub("1", "Normal", G.R.R)  ## make samples_id correct with T and N through TCGA-sampleid 
tumor.index <- grep('Tumor', G.R.R) ## get the tumor group index
tumor.index
normal.index <- grep('Normal', G.R.R)
normal.index

data.G <- rbind(G.R.R, data)  
colnames(data.G)<- paste0(c(colnames(data)), G.R.R)  ## add "T" and "N" to sample_id
data.T <- data.G[, tumor.index] 
data.N<- data.G[, normal.index]
a <- c(sample(length(data.T[1,])))
data.T.test <- data.T[, a[1: (length(data.T[1,]) %/% 2)]] ##
dim(data.T.test)
data.G.test <- cbind(data.T.test, data.N) 
dim(data.G.test)
data.T.validate <- data.T[, a[(length(data.T[1,]) %/% 2 + 1) : length(data.T[1,])]]
dim(data.T.validate)
data.G.validate <- cbind(data.T.validate, data.N)
dim(data.G.validate)


###  divided the samples into three group which are test, validate and all, so base on which group you choose
dimnames2 <- list(rownames(data.G[-1, ]), colnames(data.G[-1, ]))
data.V <- matrix(as.numeric(as.matrix(data.G[-1, ])), nrow = nrow(data.G[-1, ]), dimnames = dimnames2)
data.M <- data.V

#data.M <- data.V[apply(data.V, 1, median) > 0, ]
dim(data.M)
rns <- rownames(data.M)
cns <- colnames(data.M)
sample.id <- colnames(data.M)
rows <- nrow(data.M)
cols <- ncol(data.M)

data.L <- data.M
write.table(data.L, file = "exp_miRNA.txt", sep = "\t",row.names = TRUE, col.names = TRUE)
group <- c(colnames(data.L))
cut <- data.G[1, ]
data.L[is.na(data.L)] <- 0

###  PCA.plot
data.a <- t(as.matrix(data.L))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
pc$group <- group
pc$cut <- cut
xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file="pca_nq.tiff", width = 3600, height = 1800)
ggplot(pc, aes(PC1, PC2), main) +
  geom_point(size=15, aes(shape = cut, color = cut)) +
  scale_color_manual(values=c("green", "blue")) +
  labs(x = xlab, y = ylab, title = "LUSC") +
  theme_bw()+
  theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust = 3, face = 'bold'),
        legend.position = c(1, 1), legend.justification=c(1, 1), legend.background = element_rect(fill = 'white', colour = 'black', size = 1), 
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 3),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 3)) 

dev.off()

###  heatmap plot
### hclust
#data.L<- t(apply(data.L, 1, FUN = function(x){(x-mean(x))/sd(x)}))  # Zscore
data.G[1,][c(data.G[1,] == "Tumor")] <- 'darkblue'
data.G[1,][c(data.G[1,] == "Normal")] <- 'green'
heatmap.cols <- data.G[1,]


###  hclust

library(dendextend)
require(colorspace)

### euclidean distance, ward.D clust
tiff(file = "hclust_nq.tiff", width = 3600, height = 1800)  ###  width = 3600, height = 512
data.L.S <-t(scale(t(data.L)))
data.L.C <- hclust(d = dist(t(data.L.S )), method = 'ward.D')
data.L.D <- as.dendrogram(data.L.C)
new_order <- order.dendrogram(data.L.D)
labels(data.L.D) <- sample.id[new_order]
labels_colors(data.L.D) <- heatmap.cols[new_order]
par(cex = 1.5, font = 10, cex.axis = 5, cex.main = 5,  mar = c(3, 5, 3, 5) + 0.1)
plot(data.L.D, main = "LUSC sample id")  ###ylab = "Euclidean distance"
legend("righ", legend = c('Tumor', 'Normal'), fill = c('darkblue', 'green'), text.col = c('darkblue', 'green'), cex = 5, bty = "n")
dev.off()


###  boxplot 
gene.exp <- c(unlist(t(data.L), use.names = F))
sample.id.b <- c(rep(c(sample.id), length(data.L[,1])))
cut.b <- c(rep(c(cut), length(data.L[, 1])))
data.b <- data.frame(sample.id.b, gene.exp, cut.b, stringsAsFactors = F)
data.b$sample.id.b <- as.factor(data.b$sample.id.b)
r <- order(data.b$sample.id.b)
data.c <- data.b[r,]
DG <- c(rep("A", 96*(length(data.L[, 1]))), rep(LETTERS[2:3], each = 200*(length(data.L[, 1]))))  # base on the number of samples variated
data.d <- cbind(DG, data.c)
data.d$gene.exp <- as.numeric(c(data.d[, 3]))

tiff(file="boxplot_nq.tiff", width = 3600, height = 1800)
ggplot(data = data.d, aes(x = sample.id.b, y = gene.exp)) +
  geom_boxplot(aes(fill = cut.b), size = 1.5) +
  scale_fill_manual(values=c("green", "blue")) +
  facet_wrap(~DG, ncol = NULL, nrow = 3, scales = "free", drop = T) +
  theme_bw()+ 
  theme(axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 0.1, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 60, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 60, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust = 3, face = 'bold'),
        legend.key.size = unit(7, "cm")) 
dev.off()

#### Normalized QC

data.L.Q <- apply(data.L , 2, FUN=function(x){(x-median(x))})  ###median 0 for samples


### as above, based on which group you choose 


###  PCA.plot
data.a <- t(as.matrix(data.L.Q))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
pc$group <- group
pc$cut <- cut
xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file="pca_qc.tiff", width = 3600, height = 1800)
ggplot(pc, aes(PC1, PC3), main) +
  geom_point(size=15, aes(shape = cut, color = cut )) + 
  scale_color_manual(values=c("green", "blue")) +
  labs(x = xlab, y = ylab, title = "LUSC") +
  theme_bw()+
  theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust = 3, face = 'bold'),
        legend.position = c(1, 1), legend.justification=c(1, 1), legend.background = element_rect(fill = 'white', colour = 'black', size = 1), 
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 3),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 3))
dev.off()

###  heatmap plot
### hclust
#data.L<- t(apply(data.L, 1, FUN = function(x){(x-mean(x))/sd(x)}))  # Zscore
#data.G[1,][c(data.G[1,] == "Tumor")] <- 'red'
#data.G[1,][c(data.G[1,] == "Normal")] <- 'blue'
#heatmap.cols <- data.G[1,]


###  hclust

library(dendextend)
require(colorspace)

### euclidean distance, ward.D clust
tiff(file = "hclust_qc.tiff", width = 3600, height = 1800)  ###  width = 3600, height = 512
data.L.S <-t(scale(t(data.L.Q)))
data.L.C <- hclust(d = dist(t(data.L.S )), method = 'ward.D')
data.L.D <- as.dendrogram(data.L.C)
new_order <- order.dendrogram(data.L.D)
labels(data.L.D) <- sample.id[new_order]
labels_colors(data.L.D) <- heatmap.cols[new_order]
par(cex = 1.5, font = 10, cex.axis = 5, cex.main = 5,  mar = c(3, 5, 3, 5) + 0.1)
plot(data.L.D, main = "LUSC sample id")  ### ylab = "Euclidean distance", 
legend("righ", legend = c('Tumor','Normal'), fill = c('darkblue','green'), text.col = c('darkblue','green'), cex = 5, bty = "n")
dev.off()


###  boxplot 
gene.exp <- c(unlist(t(data.L.Q), use.names = F))
sample.id.b <- c(rep(c(sample.id), length(data.L.Q[,1])))
cut.b <- c(rep(c(cut), length(data.L.Q[, 1])))
data.b <- data.frame(sample.id.b, gene.exp, cut.b, stringsAsFactors = F)
data.b$sample.id.b <- as.factor(data.b$sample.id.b)
r <- order(data.b$sample.id.b)
data.c <- data.b[r,]
DG <- c(rep("A", 96*(length(data.L[, 1]))), rep(LETTERS[2:3], each =200*(length(data.L[, 1]))))  # base on the number of samples variated
data.d <- cbind(DG, data.c)
data.d$gene.exp <- as.numeric(c(data.d[, 3]))

tiff(file="boxplot_qc.tiff", width = 3600, height = 1800)
ggplot(data = data.d, aes(x = sample.id.b, y = gene.exp)) +
  geom_boxplot(aes(fill = cut.b), size = 1.5) +
  scale_color_manual(values=c("green", "blue")) +
  facet_wrap(~DG, ncol = NULL, nrow = 3, scales = "free", drop = T) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 0.1, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 60, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 60, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust = 3, face = 'bold'),
        legend.key.size = unit(7, "cm")) 
dev.off()

### T-test do differential expression analysis
pvalue <- function(x, y, var.equal = F, paired = F){ 
  tryCatch({
    obj <- t.test(x = x, y = y, var.equal = var.equal, paired = F)
    return(obj$p.value)
  },error=function(e){
    return(1)
  })
}



log2fpkm_based_DEGs <- function(DataSet, TumorLabel, NormalLabel, log2FC = 1, 
                                AdjustMethod= "fdr" , AdjustP = 0.05, paired = F){
  cns <- colnames(DataSet)
  TumorNames <- grep(TumorLabel, cns, value = T)
  NormalNames <- grep(NormalLabel, cns, value = T)
  mat <- DataSet[,c(sort(NormalNames), sort(TumorNames))]
  n1 <- length(NormalNames)
  n2 <- length(TumorNames)
  n <- n1+n2
  var_pvalue <- apply(mat, 1, FUN = function(x){
    ### use F test compare Tumor and Normal of mat[x,]. At last return P.value
    var.test(x[1:n1], x[(n1+1):n])$p.value 
  }) 
  na.idx <- is.na(var_pvalue)
  var_pvalue <- var_pvalue[!na.idx]
  mat <- mat[names(var_pvalue),]
  var_sam <- names(var_pvalue[var_pvalue>0.05])
  var_dif <- names(var_pvalue[var_pvalue<=0.05])
  t.test.p1 <- apply(mat[var_sam,], 1, FUN = function(x){
    pvalue(x = x[1:n1], y = x[(n1+1):n], paired = F, var.equal = T)
  })
  t.test.p2 <- apply(mat[var_dif,], 1 , FUN = function(x){
    pvalue(x = x[1:n1], y = x[(n1+1):n], paired = F, var.equal = F)
  })
  t.test.p <- c(t.test.p1, t.test.p2)
  log2fc <- apply(mat, 1, FUN=function(x){
    fc <- mean(x[(n1+1):n]) - mean(x[1:n1]); fc 
  })
  t.test.p <- t.test.p[names(log2fc)]
  t.test.p.adjut <- p.adjust(t.test.p, method = AdjustMethod)
  rns <- names(log2fc)
  diff_exp_matrix <- cbind(t.test.p, t.test.p.adjut, log2fc)
  dimnames(diff_exp_matrix) <- list(rns, c("Pvalue", toupper(AdjustMethod), "log2fc"))
  n3 <- nrow(diff_exp_matrix)
  cols <- rep(1,n3)
  de_rows <- which((diff_exp_matrix[,3] >= log2FC | diff_exp_matrix[,3] <= (-log2FC)) & diff_exp_matrix[,2] <= AdjustP)
  cols[de_rows] <- 2
  DEGs_matrix <- diff_exp_matrix[de_rows,]
  DEGs_up <- DEGs_matrix[DEGs_matrix[,3]>0,]
  DEGs_up <- DEGs_up[order(DEGs_up[,3],decreasing=T),]
  DEGs_down <- DEGs_matrix[DEGs_matrix[,3]<0,]
  DEGs_down <- DEGs_down[order(DEGs_down[,3],decreasing=F),]
  DEGs_list <- list(Matrix = diff_exp_matrix, DEGs_matrix = DEGs_matrix, 
                    DEGs_up = DEGs_up, DEGs_down = DEGs_down, colors = cols)
  return (DEGs_list)
}

log2fpkm_based_DEGs_list <- log2fpkm_based_DEGs(DataSet = data.L, TumorLabel='Tumor', NormalLabel='Normal',
                                                log2FC = 1, AdjustMethod = 'fdr', AdjustP = 0.05, paired = F)


DEG <- log2fpkm_based_DEGs_list$colors

tiff(file="Ttest_DEG_volcano.tiff", width=2500, height=4000)
qplot(x = log2fpkm_based_DEGs_list$Matrix[,3], y = -log2(log2fpkm_based_DEGs_list$Matrix[, 1]), 
      col = DEG, xlab = 'log2FC', ylab = '-log10Pvalue', size = I(15)) +
  scale_color_gradientn(colours = c("green", "blue")) +
  labs(title = "LUSC")+
  theme_bw()+
  theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), legend.key.size = unit(7, "cm"),
        axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust = 3, face = 'bold'),
        legend.background = element_rect(fill = 'white', colour = 'white', size = 1), 
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 3),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 3))
dev.off()


### classfy the differential expreesion genes which included up, down, and Z-score
Name.gene.DEG <- rownames(log2fpkm_based_DEGs_list$DEGs_matrix)
Name.gene.DEG.up <- rownames(log2fpkm_based_DEGs_list$DEGs_up)
Name.gene.DEG.down <- rownames(log2fpkm_based_DEGs_list$DEGs_down)
data.DEG <- data.L[c(match((Name.gene.DEG), rownames(data.L), nomatch = 0)), ]
#data.DEG.Z <- apply(data.DEG, 2, FUN=function(x){(x-mean(x))/sd(x)})
data.DEG.up <- data.L[c(match((Name.gene.DEG.up), rownames(data.L), nomatch = 0)), ]
data.DEG.down <- data.L[c(match((Name.gene.DEG.down), rownames(data.L), nomatch = 0)), ]



###  PCA.plot.DEG
data.a <- t(as.matrix(data.DEG))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
pc$group <- group
pc$cut <- cut
xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file="Ttest_DEG_ pca.tiff", width = 3600, height = 1800)
ggplot(pc, aes(PC1, PC2), main) +
  geom_point(size=15, aes(shape = cut, color = cut )) +
  scale_color_manual(values=c("green", "blue")) +
  #geom_text(aes(label= cut), size=0.01)+ #add text to point
  labs(x = xlab, y = ylab, title = "LUSC") +
  theme_bw()+
  theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 100, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        legend.title = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust =3, face = 'bold'),
        legend.text = element_text(colour = 'black', angle = 1, size = 100, hjust = 2, vjust = 3, face = 'bold'),
        legend.position = c(1, 1), legend.justification=c(1, 1), legend.background = element_rect(fill = 'white', colour = 'black', size = 1), 
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 3),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 3))
dev.off()

dim(log2fpkm_based_DEGs_list$DEGs_matrix)
dim(log2fpkm_based_DEGs_list$DEGs_up)
dim(log2fpkm_based_DEGs_list$DEGs_down)


write.table(data.DEG, file = 'Ttest_DEG_EXP.txt', row.names = T, col.names = T, sep = '\t', quote = F)
write.table(x = log2fpkm_based_DEGs_list$DEGs_matrix, file = 'Ttest_DEG_fp.txt', row.names = T, col.names = T, sep = '\t', quote = F)
write.table(x = log2fpkm_based_DEGs_list$DEGs_up, file = 'Ttest_DEG_up_fp.txt', row.names = T, col.names = T, sep = '\t', quote = F)
write.table(x = log2fpkm_based_DEGs_list$DEGs_down, file = 'Ttest_down_fp.txt', row.names=T, col.names=T,sep='\t',quote=F)
write.table(x = log2fpkm_based_DEGs_list$Matrix, file = 'Ttest_all_fp.txt', row.names = T, col.names = T, sep='\t',quote=F)

DEG.top <- c(row.names(log2fpkm_based_DEGs_list$DEGs_matrix[order(abs(log2fpkm_based_DEGs_list$DEGs_matrix[, 3]), decreasing = T), ])[1:100])
data.top <- data.DEG [DEG.top, ]
data.top <- apply(data.top, 2, FUN = function(x){(x-median(x))/(sd(x) + 0.0001)})

tiff(file = "heatmap_DEG.tiff", width = 3600, height = 1800)
par(cex.lab = 3, cex.axis = 2, cex.main = 3 )
heatmap.2(data.top, col = colorpanel(99, "blue", "black", "yellow"), dendrogram = "column", 
          hclustfun = function(x){hclust(x, method = 'ward.D2')}, key=T, symkey=F, key.xlab = NA,
          key.title = "Normalied expression", keysize = 1, Colv = T, Rowv = T, trace = "none", density.info = "none",
          labRow = NA, labCol = NA, cexCol = 0.5, ColSideColors = heatmap.cols, scale = "row", na.rm = T,
          margin = c(14,10), colCol = heatmap.cols, srtCol = 90) 
legend("left", legend = c('Tumor', 'Normal'), pch = 15, col = c('darkblue', 'green'), cex = 6, text.col = "black") 
dev.off() 

