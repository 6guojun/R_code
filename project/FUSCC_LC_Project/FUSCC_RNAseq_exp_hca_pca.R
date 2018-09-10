library(ggplot2)
library(tidyr)
library(dplyr)
library(gmodels)
library(plotly)
library(crosstalk)
library(sva)

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
rt_batch <- read.table(file = "batch_files.txt", header = TRUE, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
rt <- read.table(file = "FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
rt_d <- rt[-which(apply(rt[, c(405:413)], 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)

#match colnames of batch and rt 
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
rt_batch_m <- rt_batch[match(colnames(rt_l), rt_batch$samples, nomatch = NA), ] ##make sure rt_batch$samples == colnames(rt)

#remove  abnormal samples
colnames(rt_l) <- paste(matrix(unlist(strsplit(colnames(rt_l), "[_]")), ncol = 4, byrow = TRUE)[, 4])
#Num_ab <- match(c("4501T", "3126T", "4640T", "2649T", "2661T", "2967T", "1851N", "1873N"), colnames(rt_l), nomatch = NA)
#Num_ab <- match(c("4501T", "3126T", "3289T", "4640T", "3879T", "2649T", "2661T", "2967T", "1851N", "1873N"), colnames(rt_l), nomatch = NA)
#Num_ab <- match(c("4501T", "3126T", "3289T", "4640T", "3879T", "2649T", "2661T", "2967T", "1851N", "1873N", "2391N"), colnames(rt_l), nomatch = NA)
#Num_ab <- match(c("4501T", "4501N", "3126T", "3126N", "3289T", "3289N", "4640T", "4640N", "3879T", '3879N',
#                  "2649T", "2649N",  "2661T", "2661N", "2967T", "2967N", "1851N", "1851T", "1873N", "1873T", 
#                  "2391N", "2391T", "2310T", "2310N", "2254T", "2254N", "2883T", "2883N"), colnames(rt_l), nomatch = NA)
Num_ab <- 0

#add batch to colnames(rt)
colnames(rt_l) <- paste(colnames(rt_l), rt_batch_m$batch_s, sep = "_")


num_pos_A <- grep("SampleA",colnames(rt_l))
names(num_pos_A) <- rep("SampleA", each = length(num_pos_A))
num_pos_B <- grep("SampleB",colnames(rt_l))
names(num_pos_B) <- rep("SampleB", each = length(num_pos_B))

num_pos_T <- grep("T", colnames(rt_l))
names(num_pos_T) <- rep("Tumor", each = length(num_pos_T))
num_pos_N <- grep("N", colnames(rt_l))
names(num_pos_N) <- rep("Normal", each = length(num_pos_N))

rt_samp <- data.frame(cbind(colnames(rt_l), names(sort(c(num_pos_A, num_pos_B, num_pos_T, num_pos_N), decreasing = FALSE))), 
                      stringsAsFactors = FALSE)
colnames(rt_samp) <- c("samples", "group")

rt_sam_col <- rt_samp
rt_sam_col$color <- NA
rt_sam_col$color[grep("Tumor", rt_samp$group)] <- "red"
rt_sam_col$color[grep("Normal", rt_samp$group)] <- "blue"
rt_sam_col$color[grep("SampleA", rt_samp$group)] <- "purple"
rt_sam_col$color[grep("SampleB", rt_samp$group)] <- "green"
rt_sam_col$batch_s <- rt_batch_m$batch_s
colnames(rt_sam_col) <- c("samples", "group", "color", "batch_s")


#################################################exp_hca_pca############################################


### make a boxplot of all samples form gene level
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
MakBoxPlot(rt_l, rt_sam_col, theme_E, type = "tiff", width = 2048, height = 1024)


###sampleA and B were removed 
#rt_l_d <- rt_l[, -c(405:413, Num_ab)]
#rt_sam_col_d <- rt_sam_col[-c(405:413, Num_ab), ]
rt_l_d <- rt_l[, c(405:413)]
rt_sam_col_d <- rt_sam_col[c(405:413), ]
#rt_l_d <- rt_l[, Num_ab]
#rt_sam_col_d <- rt_sam_col[Num_ab, ]

### do hca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_l_d, rt_sam_col = rt_sam_col_d, "LC_batch1", type = "tiff", method = 'ward.D', cex = 1, width = 3600, height = 1024)

### make a pca plot of all samples 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l_d, rt_sam = rt_sam_col_d, nam = "ALL", ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)

data.a <- t(as.matrix(rt_l_d))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_sam_col_d$samples, nomatch = NA)
pc$group <-  rt_sam_col_d$group[num_pos]
pc$color <-  rt_sam_col_d$color[num_pos]

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "ALL_PCA_text.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group)) +  scale_colour_manual(values = unique(pc$color)) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)


###combat remove batch effected
pheno <- rt_sam_col_d[, c(2, 4)]
row.names(pheno) <- rt_sam_col_d$samples
batch <- pheno$batch_s
modcombat <- model.matrix(~1, data = pheno)
rt_cBat <- ComBat(dat = rt_l_d, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

### make a boxplot of all samples form gene level
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
MakBoxPlot(rt_cBat, rt_sam_col_d, theme_E, type = "tiff", width = 3600, height = 1600)

### do hca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_cBat, rt_sam_col = rt_sam_col_d, "LC_ALL", type = "tiff", method = 'complete', cex = 1.8, width = 3600, height = 1600)


### make a pca plot of all samples 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_cBat, rt_sam = rt_sam_col_d, nam = "ALL", ThemePCA = theme_A, type = "tiff", width = 1024, height = 512)

data.a <- t(as.matrix(rt_cBat))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_sam_col_d$samples, nomatch = NA)
pc$group <-  rt_sam_col_d$group[num_pos]
pc$color <-  rt_sam_col_d$color[num_pos]

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "FUSCC_LC_PCA_text_batch1.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group)) +  scale_colour_manual(values = unique(pc$color)) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)

##########################################DEG analysis######################################################

###DEG
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/DEG_T_test.R")
DEGPvalFC <- DEGTest("Tumor", "Normal", rt_l_d, rt_sam_col_d, var.equal = FALSE, paired = FALSE, meas = "DEGPvalFC")
DEG_S <- DEGPvalFC[order(DEGPvalFC$log2FC), ]
rt_l_d_d <- rt_l_d[row.names(DEG_S), ]

###top100 gene list
#gene_list <- c(row.names(DEG_S)[1 :50], row.names(DEG_S)[5354: 5854])
#rt_l_d_t <- rt_l_d[gene_list, ]

###do heatmap
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/heatmap.R")
MakHeatmap(rt_l_d_d, rt_sam = rt_sam_col_d, "Tumor", "Normal", type = "tiff", width = 1500, height = 1024)


###make a HCA
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_l_d_d, rt_sam_col = rt_sam_col_d, "DEG", type = "tiff", method = 'ward.D2', cex = 1, width = 3600, height = 1024)


### PCA analysis
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l_d_d, rt_sam = rt_sam_col_d, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)


data.a <- t(as.matrix(rt_l_d_d))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_sam_col_d$samples, nomatch = NA)
pc$group <-  rt_sam_col_d$group[num_pos]
pc$color <-  rt_sam_col_d$color[num_pos]

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "DEG_PCA_text.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group)) +  scale_colour_manual(values = unique(pc$color)) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)


########################################DEG_test_remove#####################################
###DEG gene list
gene_list <- row.names(DEG_S)
rt_l_d_t <- rt_l_d[gene_list, ]


###do heatmap
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/heatmap.R")
MakHeatmap(rt_l_d_t, rt_sam = rt_sam_col_d, "Tumor", "Normal", type = "tiff", width = 1500, height = 1024)


###make a HCA
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_l_d_t, rt_sam_col = rt_sam_col_d, "DEG", type = "tiff", method = 'ward.D2', cex = 1, width = 3600, height = 1024)


### PCA analysis
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l_d_t, rt_sam = rt_sam_col_d, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)


data.a <- t(as.matrix(rt_l_d_t))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_sam_col_d$samples, nomatch = NA)
pc$group <-  rt_sam_col_d$group[num_pos]
pc$color <-  rt_sam_col_d$color[num_pos]

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "DEG_PCA_text.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group)) +  scale_colour_manual(values = unique(pc$color)) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)



#########################################gene id convert############################################
row.names(rt_l_d) <-  matrix(unlist(strsplit(row.names(rt_l_d), "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg <- bitr(row.names(rt_l_d), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
rt_l_d_s <- cbind(row.names(rt_l_d), rt_l_d[, -397])
rt_l_d_c <- merge(eg, rt_l_d_s, by = 1)
rt_l_d_c_d <- rt_l_d_c[, -1]
write.table(rt_l_d_c_d, file = "FUSCC_LC_RNAseq_FPKM_log_symbel_R.txt", row.names = FALSE, col.names = TRUE, sep = "\t")




#########################################DEG_more_analysis#########################################
setwd("/Volumes/Macintosh HD/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp/DEG_remove8")
rt_DEG <- read.table(file = "Tumor_Normal_DEG_Exp.txt", header = TRUE, sep = "\t", row.names  = 1)


###make all boxplot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/Mutiple_violin_boxplot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
lapply(row.names(DEGExp), MulVBplot, DEGExp, rt_sam = NULL, GP1 = "WT", GP2 = "Mut", theme_D, var.equal = FALSE)