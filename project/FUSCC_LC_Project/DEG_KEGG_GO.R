library(ggplot2)
library(tidyr)
library(gmodels)
#library(crosstalk)
#library(sva)

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")

rt_batch <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp/batch_files.txt", header = TRUE, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
rt <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp/FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
rt_d <- rt[-which(apply(rt[, -c(405:413)], 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp/ALL_DEG")
#match colnames of batch and rt 
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
rt_batch_m <- rt_batch[match(colnames(rt_l), rt_batch$samples, nomatch = NA), ] ##make sure rt_batch$samples == colnames(rt)

#remove  abnormal samples
colnames(rt_l) <- paste(matrix(unlist(strsplit(colnames(rt_l), "[_]")), ncol = 4, byrow = TRUE)[, 4])
#Num_ab <- match(c("4501T", "3126T", "4640T", "2649T", "2661T", "2967T", "1851N", "1873N"), colnames(rt_l), nomatch = NA)
#Num_ab <- match(c("4501T", "3126T", "3289T", "4640T", "3879T", "2649T", "2661T", "2967T", "1851N", "1873N"), colnames(rt_l), nomatch = NA)
#Num_ab <- match(c("4501T", "3126T", "3289T", "4640T", "3879T", "2649T", "2661T", "2967T", "1851N", "1873N", "2391N"), colnames(rt_l), nomatch = NA)
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

###sampleA and B were removed 
rt_l_d <- rt_l[, -c(405:413, Num_ab)]
rt_sam_col_d <- rt_sam_col[-c(405:413, Num_ab), ]

##########################################DEG analysis######################################################

###DEG
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/DEG_T_test.R")
DEGPvalFC <- DEGTest("Tumor", "Normal", rt_l_d, rt_sam_col_d, var.equal = FALSE, paired = FALSE, meas = "DEGPvalFC")
DEGExp <- DEGTest("Tumor", "Normal", rt_l_d, rt_sam_col_d, var.equal = FALSE, paired = FALSE, meas = "DEGExp")
AllPvalFC <- DEGTest("Tumor", "Normal", rt_l_d, rt_sam_col_d, var.equal = FALSE, paired = FALSE, meas = "AllPvalFC")
#DEGTest("Tumor", "Normal", rt_l_d, rt_sam_col_d, var.equal = FALSE, paired = FALSE, meas = "")
DEG_S <- DEGPvalFC[order(DEGPvalFC$log2FC), ]
rt_l_d_d <- rt_l_d[row.names(DEG_S), ]

###top100 gene list
gene_list <- c(row.names(DEG_S)[1 :50], row.names(DEG_S)[(length(DEG_S$Pval) - 50) : length(DEG_S$Pval)])
rt_l_d_t <- rt_l_d[gene_list, ]

###do heatmap
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/heatmap.R")
MakHeatmap(rt_l_d_t, rt_sam = rt_sam_col_d, "Tumor", "Normal", type = "tiff", width = 1500, height = 1024)

### volcano
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/volcano_plot.R")
MakVolPlot(AllPvalFC, "Tumor", "Normal", 2, theme_D)

###boxplot
setwd("/mnt/home_91TB/home/shangjun/project_level4/FUSCC_LC_RNAseq_ALL_Hisat2_0_5/ALL_DEG/boxplot_DEG")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script//Mutiple_violin_boxplot.R")
lapply(row.names(DEGExp), MulVBplot, DEGExp, rt_sam = NULL, GP1 = "Tumor", GP2 = "Normal", theme_D, var.equal = FALSE)

###KEGG and GO
setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp/ALL_KEGG_GO")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/KEGG_GO.R")
AllPvalFC$gnam <- matrix(unlist(strsplit(rownames(AllPvalFC), "[.]")), byrow = TRUE, ncol = 2)[, 1]
MakKEGO(AllPvalFC, ont = "MF", data_type = "gene_list", Gene_ID = "ENSEMBL", fc_v = 1, 
        OrgDb = "org.Hs.eg.db", organism = "hsa")
                                                            
