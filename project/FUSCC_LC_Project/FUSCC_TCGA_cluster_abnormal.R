library(ggplot2)
library(tidyr)
library(dplyr)
library(gmodels)
library(plotly)
library(crosstalk)
library(sva)


setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
rt_all <- read.table(file = "FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
colnames(rt_all) <- matrix(unlist(strsplit(colnames(rt_all), "[.]")), ncol = 2, byrow = TRUE)[, 2]
colnames(rt_all) <- matrix(unlist(strsplit(colnames(rt_all), "[_]")), ncol = 4, byrow = TRUE)[, 4]

rt_e <- rt_all[, c("3126T", "4941T", "2649T", "4501T", "1873N", "4402T")]

num_pos_T <- grep("T", colnames(rt_e))
names(num_pos_T) <- rep("FUSCC_T", each = length(num_pos_T))
num_pos_N <- grep("N", colnames(rt_e))
names(num_pos_N) <- rep("FUSCC_N", each = length(num_pos_N))

rt_samp_F <- data.frame(cbind(colnames(rt_e), names(sort(c(num_pos_T, num_pos_N), decreasing = FALSE))), 
                        stringsAsFactors = FALSE)
colnames(rt_samp_F) <- c("samples", "group")


###get TCGA data
rt_tcga <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/TCGA_LUAD_FPKM.txt", header = TRUE, row.names = 1, sep = "\t", 
                      stringsAsFactors = FALSE)

sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "[.]")), ncol = 7, byrow = TRUE)[, 4]

tcga_pos_T1 <- grep("01", sam_G)
names(tcga_pos_T1) <- rep("TCGA_T", each = length(tcga_pos_T1))
tcga_pos_T2 <- grep("02", sam_G)
names(tcga_pos_T2) <- rep("TCGA_T", each = length(tcga_pos_T2))
tcga_pos_N <- grep("11", sam_G)
names(tcga_pos_N) <- rep("TCGA_N", each = length(tcga_pos_N))
rt_samp_T <- data.frame(cbind(colnames(rt_tcga), names(sort(c(tcga_pos_T1, tcga_pos_T2, tcga_pos_N), decreasing = FALSE))), 
                        stringsAsFactors = FALSE)
colnames(rt_samp_T) <- c("samples", "group")


###combine FUSCC and TCGA
rt_e <- rt_e[match(row.names(rt_tcga), row.names(rt_e), nomatch = NA), ]
rt_all_tcga <- cbind(rt_e, rt_tcga)
rt_d <- rt_all_tcga[-which(apply(rt_all_tcga, 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
colnames(rt_l) <- colnames(rt_d)

rt_sam <- rbind(rt_samp_F, rt_samp_T)
rt_sam$batch <- c(rep(1, length(grep("FUSCC", rt_sam$group))), rep(2, length(grep("TCGA", rt_sam$group))))

rt_sam_col <- rt_sam
rt_sam_col$color <- NA
rt_sam_col$color[grep("FUSCC_T", rt_sam$group)] <- "red"
rt_sam_col$color[grep("FUSCC_N", rt_sam$group)] <- "blue"
rt_sam_col$color[grep("TCGA_T", rt_sam$group)] <- "purple"
rt_sam_col$color[grep("TCGA_N", rt_sam$group)] <- "green"

###ComBat remove batch which come from TCGA and FUSCC
pheno <- rt_sam_col[, c('group', "batch")]
row.names(pheno) <- rt_sam_col$samples
batch <- pheno$batch
modcombat <- model.matrix(~1, data = pheno)
rt_l_combat = ComBat(dat = rt_l, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

### do hca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_l_combat, rt_sam_col = rt_sam_col, "LC_batch1", type = "tiff", method = 'ward.D', cex = 1.8, width = 1024, height = 716)

### make a pca plot of all samples 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l, rt_sam = rt_sam_col, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)
MakPCA(rt_l_combat, rt_sam = rt_sam_col, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)


data.a <- t(as.matrix(rt_l))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_sam$samples, nomatch = NA)
pc$group <-  rt_sam$group

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "FUSCC_TCGA_LUAD_PCA_ALL.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group), main) + 
  geom_point(size = 5) + #geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)