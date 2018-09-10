library(survival)
library(ggplot2)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(annotables)
library(pathview)
library(RColorBrewer)
library(scales)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/clinical_analysis_results")

cancer="GBM"

####################################################
#get cancer exprssion 
#get the number of mutation and neo-ag
####################################################

####cancer expression and convert
rt <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                 sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt$Ensembl_ID <- matrix(unlist(strsplit(rt$Ensembl_ID, "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(rt$Ensembl_ID, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
rt_sym <- merge(eg, rt, by = 1)
row.names(rt_sym) <- paste(rt_sym$ENSEMBL, rt_sym$SYMBOL, sep = "_")
rt_sym <- rt_sym[, -c(1, 2)]

#for(i in 1:61146){
#  logical_test <- all(rt_sym[grep(rt$Ensembl_ID[i], rt_sym$ENSEMBL), ][, -c(1, 2)] == rt[i, -1])
#  print(logical_test)
#  print(i)
#}



source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")
rt_T <- split_tcga_tn(rt_sym, sam_type = "tumor", split_type = "[.]")
#rt_N <- split_tcga_tn(rt_sym[, -1], sam_type = "normal")

#============================================================
###CTLA4 entrezid is 1493
#CTLA4 <- as.numeric(rt_T[which(rt_sym$SYMBOL == "CTLA4"), ])
#names(CTLA4) <- colnames(rt_T)
#29126
#PDL1 <- as.numeric(rt_T[which(rt_sym$SYMBOL == "CD274"), ])
#names(PDL1) <- colnames(rt_T)
#5133
#PD1  <- as.numeric(rt_T[which(rt_sym$SYMBOL == "PDCD1"), ])
#names(PD1) <- colnames(rt_T)
#50943
#FOXP3 <- as.numeric(rt_T[which(rt_sym$SYMBOL == "FOXP3"), ])
#names(FOXP3) <- colnames(rt_T)
#============================================================

###clinical ddta
rt_cli <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/phenotype/TCGA-", cancer, ".GDC_phenotype.tsv", sep = ""),
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE, fill = TRUE, quote = "", na.strings = "NA")
sam_cli <- matrix(unlist(strsplit(rt_cli$submitter_id.samples, "-")), ncol = 4, byrow = TRUE)
rt_cli$submitter_id.samples <- paste(sam_cli[, 1], sam_cli[, 2], sam_cli[, 3], sep = "-")
pari_sam <- match(colnames(rt_T), rt_cli$submitter_id.samples, nomatch = 0)
rt_cli_M <- rt_cli[pari_sam, ]
rt_cli_age <- as.numeric(rt_cli_M$age_at_initial_pathologic_diagnosis) 
rt_cli_age[which(rt_cli_age == "")] = NA
rt_cli_gender <- rt_cli_M$gender.demographic
rt_cli_gender[which(rt_cli_gender == "")] = NA
rt_cli_stage <- rt_cli_M$tumor_stage.diagnoses
rt_cli_stage[which(rt_cli_stage == "")] = NA
rt_cli_race <- rt_cli_M$race.demographic
rt_cli_race[which(rt_cli_race == "")] = NA

###make a boxplot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
sam_gender <- data.frame(cbind(rt_cli_M$submitter_id.samples, rt_cli_gender), stringsAsFactors = FALSE)
colnames(sam_gender) <- c("samples", "group")

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/clinical_analysis_results/gender")
for(nam in row.names(rt_T)){
  rt_gene <- rt_T[nam, ]
  rt_box_geneder <- conver_exp_data(rt_gene, sam_gender, col_names = c("samples_id", "group", "gene_id", nam))
  MakBoxPlot(nam, "group", rt_box_geneder, width = 6, height = 6, theme_B) 
}

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/clinical_analysis_results/race")
sam_race <- data.frame(cbind(rt_cli_M$submitter_id.samples, rt_cli_race), stringsAsFactors = FALSE)
colnames(sam_race) <- c("samples", "group")
for(nam in row.names(rt_T)){
  rt_gene <- rt_T[nam, ]
  rt_box_race <- conver_exp_data(rt_gene, sam_race, col_names = c("samples_id", "group", "gene_id", nam))
  MakBoxPlot(nam, "group", rt_box_race, width = 8, height = 10, theme_E) 
}

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/clinical_analysis_results/stage")
sam_stage <- data.frame(cbind(rt_cli_M$submitter_id.samples, rt_cli_stage), stringsAsFactors = FALSE)
colnames(sam_stage) <- c("samples", "group")
for(nam in row.names(rt_T)){
  rt_gene <- rt_T[nam, ]
  rt_box_stage <- conver_exp_data(rt_gene, sam_stage, col_names = c("samples_id", "group", "gene_id", nam))
  MakBoxPlot(nam, "group", rt_box_stage, width = 6, height = 6, theme_B) 
}



