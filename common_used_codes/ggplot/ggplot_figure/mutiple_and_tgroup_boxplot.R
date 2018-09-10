library(ggplot2)
library(clusterProfiler)
library(circlize)
library(org.Hs.eg.db)
library(pathview)

color_type_table <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/color.txt", sep = "\t", stringsAsFactors = FALSE)
color_type <- color_type_table$V2
  


setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/gene_all_cancer_exp")

cancer="GBM"
cancer_list <- c("GBM", "LUAD")
###read all tcga cancer expressino table and merge
###cancer list



source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")

GetBigMat <- function(cancer_list){
  cancer <- cancer_list[1]
  rt_on <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""), 
                      sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  rt_T <- split_tcga_tn(rt_on, sam_type = "tumor", split_type = "[.]")
  rt_N <- split_tcga_tn(rt_on, sam_type = "normal", split_type = "[.]")
  rt_on_T_m <- data.frame(cbind(cancer, t(rt_T)), stringsAsFactors = FALSE)
  rt_on_N_m <- data.frame(cbind(cancer, t(rt_N)), stringsAsFactors = FALSE)
  
  for(i in 2: length(cancer_list)){
    ###add the other cancer types
    cancer_add <- cancer_list[i]
    rt_add <-read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer_add, "/RNA_seq/TCGA-", cancer_add, ".htseq_fpkm.tsv", sep = ""), 
                        sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rt_add_T <- split_tcga_tn(rt_add, sam_type = "tumor", split_type = "[.]")
    rt_add_N <- split_tcga_tn(rt_add, sam_type = "normal", split_type = "[.]")
    rt_add_T_m <- data.frame(cbind(cancer_add, t(rt_add_T)), stringsAsFactors = FALSE)
    rt_add_N_m <- data.frame(cbind(cancer_add, t(rt_add_N)), stringsAsFactors = FALSE)
    colnames(rt_add_T_m)[1] <- "cancer"
    colnames(rt_add_N_m)[1] <- "cancer"
    if(all(colnames(rt_add_T_m) == colnames(rt_add_T_m))){
      rt_on_T_m <- data.frame(rbind(rt_on_T_m, rt_add_T_m), stringsAsFactors = FALSE)
      rt_on_N_m <- data.frame(rbind(rt_on_N_m, rt_add_N_m), stringsAsFactors = FALSE)
    } else {
      stop('the gene names and order of all cancer must be same')
    }
  }
 return(list(rt_on_T_m, rt_on_N_m))
}

rt_all_list <- GetBigMat(cancer_list)
rt_cancer_all <- rt_all_list[[1]]
rt_cancer_all$group <- "tumor"
rt_normal_all <- rt_all_list[[2]]
rt_normal_all$group <- "normal"
rt_all <- rbind(rt_cancer_all, rt_normal_all)

ensemble_gene <- matrix(unlist(strsplit(colnames(rt_all)[-c(1, 60485)], "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(ensemble_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
symbol_gene <- eg$SYMBOL[match(ensemble_gene, eg$ENSEMBL, nomatch = 0)]
colnames(rt_all) <- c("cancer_type", paste(ensemble_gene, symbol_gene, sep = "_"), "group")


###make a boxlot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
for(nam in colnames(rt_all)[-1]){
  rt_box <- rt_all[, c("cancer_type", "group", nam)]
  rt_box[, nam] <- as.numeric(rt_box[, nam])
  pdf(file = paste(nam, "_all_cancer_exp.pdf", sep = ""), 25, 10)
  p <- ggplot(rt_box, aes(x = reorder(cancer_type, rt_box[, nam], FUN=median), y = rt_box[, nam], fill = group)) +  
    geom_boxplot() + scale_fill_manual(values= color_type[c(25, 33)]) + ylab("log2(expression + 1)") + theme_E
  print(p)
  dev.off()
}


CountPal <- function(cnam, rt_all){
  rt_cnam <-rt_all[grep(cnam, rt_all[, "cancer_type"]), ]
  nam_one <- colnames(rt_all)[-c(1, 60485)][1]
  gene_exp_on <- rt_cnam[, c("cancer_type", "group", nam_one)]
  T_exp_one <- as.numeric(gene_exp_on[, nam_one][grep("tumor",rt_cnam[, "group"])])
  N_exp_one <- as.numeric(gene_exp_on[, nam_one][grep("normal",rt_cnam[, "group"])])
  pval_on_list <-  t.test(N_exp_one, T_exp_one)
  pval_on <- c(cnam, nam_one, pval_on_list$p.value)
  
  for(gnam in colnames(rt_all)[-c(1, 2, 60485)][1:10]){
    gene_exp_add <- rt_cnam[, c("cancer_type", "group", gnam)]
    T_exp_add <- as.numeric(gene_exp_add[, gnam][grep("tumor",rt_cnam[, "group"])])
    N_exp_add <- as.numeric(gene_exp_add[, gnam][grep("normal",rt_cnam[, "group"])])
    pval_add_list <- t.test(N_exp_add, T_exp_add)
    pval_add <- c(cnam, gnam, pval_add_list$p.value)
    pval_on <- rbind(pval_on, pval_add)
  }
  return(pval_on)
}

pval_list <- lapply(cnam, CountPal, rt_all)
pval_table <- do.call(cbind, pval_list)
write.table(pval_table, file = "all_cancer_TN_pvalue.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
