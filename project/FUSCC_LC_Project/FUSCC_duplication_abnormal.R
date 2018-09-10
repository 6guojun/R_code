library(ggplot2)
library(tidyr)
library(dplyr)
library(gmodels)
library(plotly)
library(crosstalk)
library(sva)

setwd("/Users/stead/Desktop/FUSCC_LC_QC")
rt_maqc <- read.table(file = "LC_MAQC_FPKM.txt", header = TRUE, sep = "\t")
rt_d <- rt_maqc[-which(apply(rt_maqc, 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
colnamesrt_l <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
row.names(rt_l) <-  matrix(unlist(strsplit(row.names(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg <- bitr(row.names(rt_l), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
rt_l_d_s <- cbind(row.names(rt_l), rt_l)
rt_l_d_c <- merge(eg, rt_l_d_s, by = 1)
rt_l_d_c_d <- rt_l_d_c[, -c(1:2)]

rt_taq <- read.table(file = "/Users/stead/Desktop/RNA_SeqTest/Taqman_raw_data/MAQC_TAQ_16r1044_ANA_gene_20150804.txt", header = TRUE, sep = "\t")
num_taq <- match(rt_taq$GeneName, rt_l_d_c$SYMBOL, nomatch = 0)
rt_l_taq <- rt_l_d_c_d[num_taq, ]

###make a correlation among 9 MAQC samples
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/Mut_correlation.R")
MakCorPlot(rt_l_taq, Tnam = "MAQC_LC_A_B", height = 1024, width = 1024, cex.axis = 2)



#######################################duplication abnoraml samples#####################################

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
rt_batch <- read.table(file = "batch_files.txt", header = TRUE, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
rt <- read.table(file = "FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
rt_d <- rt[-which(apply(rt[, c(405:413)], 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)

#match colnames of batch and rt 
colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
rt_batch_m <- rt_batch[match(colnames(rt_l), rt_batch$samples, nomatch = NA), ] ##make sure rt_batch$samples == colnames(rt)

#remove  high duplication samples
colnames(rt_l) <- paste(matrix(unlist(strsplit(colnames(rt_l), "[_]")), ncol = 4, byrow = TRUE)[, 4])
Num_ab <- match(c("3412N", "4282T", "1229N", "1906N", "1553N", "2801N", "2383N", "3157N", "4501N", 
                  "4501T", "3040N", "4553N", "1653N", "1850N", "2373N", "3146N", "3146T", "3321T", 
                  "2887N", "3218N", "3466T", "3875T", "4428N", "1964T", "2569N", "2298N", "2364T", 
                  "2926T", "3347T", "3426N", "4402N", "4402T", "4597N", "4597T", "2254N", "2373T", 
                  "2753N", "1877T", "3899N", "3899T", "3303T"), colnames(rt_l), nomatch = NA)


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
rt_l_d <- rt_l[, c(Num_ab, 405:413)]
rt_sam_col_d <- rt_sam_col[c(Num_ab, 405:413), ]

rt_l_d_n <- rt_l[, c(2:25, 112: 127, 140: 149)]
rt_sam_col_n <- rt_sam_col[c(2:25, 112: 127, 140: 149), ]

### make a boxplot of all samples form gene level
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
MakBoxPlot(rt_l_d_n, rt_sam_col_n, theme_E, type = "tiff", width = 2048, height = 1024)


###sampleA and B were removed 
#rt_l_d <- rt_l[, -c(405:413, Num_ab)]
#rt_sam_col_d <- rt_sam_col[-c(405:413, Num_ab), ]

#rt_l_d <- rt_l[, Num_ab]
#rt_sam_col_d <- rt_sam_col[Num_ab, ]
