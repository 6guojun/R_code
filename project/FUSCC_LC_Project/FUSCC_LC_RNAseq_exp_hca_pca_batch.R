###FUSCC_LC RNA-seq data
###exp, hca, pca
###base on each batch



library(ggplot2)
library(tidyr)
library(dplyr)
library(gmodels)
library(plotly)
library(crosstalk)

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
rt_batch <- read.table(file = "batch_files.txt", header = TRUE, row.names = NULL, sep = "\t", quote = "")
rt <- read.table(file = "FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
rt_d <- rt[-which(apply(rt[, -1], 1, mean) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[_]")), ncol = 4, byrow = TRUE)[, 4]

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
colnames(rt_sam_col) <- c("samples", "group", "color")



### make a boxplot of all samples form gene level
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
MakBoxPlot(rt_l, rt_sam_col, theme_E, type = "tiff", width = 1600, height = 1024)


rt_l_d <- rt_l[, -1]
rt_sam_col_d <- rt_sam_col[-1, ]

### do hca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_l_d, rt_sam_col = rt_sam_col_d, "LC_batch1", type = "tiff", method = 'ward.D', cex = 1.8, width = 1024, height = 716)

### make a pca plot of all samples 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l_d, rt_sam = rt_sam_col_d, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)

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

tiff(file = "FUSCC_LC_PCA_text_batch1.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group)) +  scale_colour_manual(values = unique(pc$color)) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)


setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
rt_batch1 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch1/exp/FUSCC_LC_FPKM_batch1.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch2 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch2/exp/FUSCC_LC_FPKM_batch2.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch3 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch3/exp/FUSCC_LC_FPKM_batch3.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch4 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch4/exp/FUSCC_LC_FPKM_batch4.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch5 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch5/exp/FUSCC_LC_FPKM_batch5.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch6 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch6/exp/FUSCC_LC_FPKM_batch6.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch7 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch7/exp/FUSCC_LC_FPKM_batch7.txt", header = TRUE, row.names = 1, sep = "\t")
rt_batch8 <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch8/exp/FUSCC_LC_FPKM_batch8.txt", header = TRUE, row.names = 1, sep = "\t")

rt_all <- data.frame(cbind(rt_batch1, rt_batch2, rt_batch3, rt_batch4, rt_batch5, rt_batch6, 
                           rt_batch7, rt_batch8), stringsAsFactors = FALSE)

num_pos_A <- grep("SampleA",colnames(rt_all))
names(num_pos_A) <- rep("SampleA", each = length(num_pos_A))

rt_d <- rt_all[-which(apply(rt_all[, -num_pos_A], 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[_]")), ncol = 4, byrow = TRUE)[, 4]


num_pos_T <- grep("T", colnames(rt_l))
names(num_pos_T) <- rep("Tumor", each = length(num_pos_T))
num_pos_N <- grep("N", colnames(rt_l))
names(num_pos_N) <- rep("Normal", each = length(num_pos_N))

rt_samp <- data.frame(cbind(colnames(rt_l), names(sort(c(num_pos_A, num_pos_T, num_pos_N), decreasing = FALSE))), 
                      stringsAsFactors = FALSE)
colnames(rt_samp) <- c("samples", "group")

rt_sam_col <- rt_samp
rt_sam_col$color <- NA
rt_sam_col$color[grep("Tumor", rt_samp$group)] <- "red"
rt_sam_col$color[grep("Normal", rt_samp$group)] <- "blue"
rt_sam_col$color[grep("SampleA", rt_samp$group)] <- "cyan"
colnames(rt_sam_col) <- c("samples", "group", "color")

### make a boxplot of all samples form gene level
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
MakBoxPlot(rt_l, rt_sam_col, theme_E, type = "tiff", width = 2600, height = 1600)

### do hca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
rt_l_d <- rt_l[, -num_pos_A]
rt_sam_col <- rt_samp
rt_sam_col$color <- NA
rt_sam_col$color[grep("Tumor", rt_samp$group)] <- "blue"
rt_sam_col$color[grep("Normal", rt_samp$group)] <- "red"
rt_sam_col$color[grep("SampleA", rt_samp$group)] <- "cyan"
rt_sam_col <- rt_sam_col[-num_pos_A, ]
MakHCA(rt_l_d, rt_sam_col = rt_sam_col, "LC_batch1", type = "tiff", method = 'ward.D', cex = 1.8, width = 2600, height = 1600)

### make a pca plot of all samples 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l_d, rt_sam = rt_sam_col, ThemePCA = theme_D, type = "tiff", width = 2048, height = 1024)

data.a <- t(as.matrix(rt_l))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_samp$samples, nomatch = NA)
pc$group <-  rt_samp$group[num_pos]

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "FUSCC_LC_PCA_text_batch1.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group), main) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)