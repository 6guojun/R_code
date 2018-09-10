library(clusterProfiler)
library(ggplot2)
setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/correlation_and_go_cluster")

cancer="GBM"

####################################################
#get cancer exprssion 
#get the number of mutation and neo-ag
####################################################

####cancer expression and convert
rt <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                 sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
row.names(rt) <- matrix(unlist(strsplit(row.names(rt), "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(row.names(rt), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
rt_sym <- merge(eg, cbind(row.names(rt), rt), by = 1)[, -1] 
###there are many duplication genes in rt_sym

source("/Users/stead/Desktop/PD-L1_and_TMI_type/get_tcga_sample_type.R")
rt_T <- split_tcga_tn(rt_sym[, -1], sam_type = "tumor", split_type = "[.]")
#rt_N <- split_tcga_tn(rt_sym[, -1], sam_type = "normal")

CTLA4 <- as.numeric(rt_T[which(rt_sym$SYMBOL == "CTLA4"), ])
names(CTLA4) <- colnames(rt_T)

Ag_pre <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/antigen_presentation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cell_ad <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/cell_adhesion.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
co_in <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/co_inhibitor.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
co_st <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/co_stimular.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ligand <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/Ligand.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
others <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/others.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
receptor <- read.table(file= "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_checkpoint/receptor.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


Ag_pre_exp <- rt_T[match(c("CTLA4", Ag_pre$Antigen_presentation), rt_sym$SYMBOL, nomatch = NA), ]
row.names(Ag_pre_exp) <- c("CTLA4", Ag_pre$Antigen_presentation)

cell_ad_exp <- rt_T[match(c("CTLA4", cell_ad$cell.adhesion), rt_sym$SYMBOL, nomatch = NA), ]
row.names(cell_ad_exp) <- c("CTLA4", cell_ad$cell.adhesion)

co_in_exp <- rt_T[match(c("CTLA4", co_in$co.inhibitor), rt_sym$SYMBOL, nomatch = NA), ]
row.names(co_in_exp) <- c("CTLA4", co_in$co.inhibitor)

co_st_exp <- rt_T[match(c("CTLA4", co_st$co.stimular), rt_sym$SYMBOL, nomatch = NA), ]
row.names(co_st_exp) <- c("CTLA4", co_st$co.stimular)

ligand_exp <- rt_T[match(c("CTLA4", ligand$Ligand), rt_sym$SYMBOL, nomatch = NA), ]
row.names(ligand_exp) <- c("CTLA4", ligand$Ligand)

others_exp <- rt_T[match(c("CTLA4", others$others), rt_sym$SYMBOL, nomatch = NA), ]
row.names(others_exp) <- c("CTLA4", others$others)

receptor_exp <- rt_T[match(receptor$receptor, rt_sym$SYMBOL, nomatch = NA), ]
row.names(receptor_exp) <- receptor$receptor


source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/theme_cor_heatmap.R") 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/correlation_heatmap.R")
cor_Ag_pre <- t(Ag_pre_exp)
cell_ad_exp <- t(cell_ad_exp)
co_in_exp <- t(co_in_exp)
co_st_exp <- t(co_st_exp)
ligand_exp <- t(ligand_exp)
others_exp <- t(others_exp)
receptor_exp <- t(receptor_exp)

mak_cor_heatmap(cor_Ag_pre, "cor_Ag_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)
mak_cor_heatmap(cell_ad_exp, "cell_ad_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)
mak_cor_heatmap(co_in_exp, "co_in_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)
mak_cor_heatmap(co_st_exp, "Aco_st_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)
mak_cor_heatmap(ligand_exp, "ligand_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)
mak_cor_heatmap(others_exp, "others_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)
mak_cor_heatmap(receptor_exp, "receptor_cor_heatmap", width = 10, height = 10, cor_text_size = 4, theme_cor_heatmap)




