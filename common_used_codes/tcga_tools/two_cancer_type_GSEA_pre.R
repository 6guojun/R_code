setwd('/Users/stead/Desktop/GSEA/two_cancer_type')

cancer1 <- 'GBM'
cancer2 <- 'LGG'


rt_c1 <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer1, "/RNA_seq/TCGA-", cancer1, ".htseq_fpkm.tsv", sep = ""), 
                    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_c2 <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer2, "/RNA_seq/TCGA-", cancer2, ".htseq_fpkm.tsv", sep = ""), 
                    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)


###get tumor samples
source("/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/get_tcga_sample_type.R")
rt_c1_T <- split_tcga_tn(rt_c1, sam_type = "tumor", split_type = "[.]")
rt_c2_T <- split_tcga_tn(rt_c2, sam_type = "tumor", split_type = "[.]")
  
###merge GBM and LGG
rt_M <- cbind(rt_c1_T, rt_c2_T)

###convert the expression table with ensemble id to gene symbol
source('/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/Ensemble2SymbelM.R')
rt_M_S <- E2SM(rt_M)
rt_M_S_d <- rt_M_S[apply(rt_M_S, 1, function(x){median(x) > 0}), ]
rt_M_S_l <- data.frame(apply(rt_M_S_d, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
rt_M_S_l <- cbind(row.names(rt_M_S_l), rt_M_S_l)
colnames(rt_M_S_l)[1] <- 'gene_id'
write.table(rt_M_S_l, file = paste(cancer1, cancer2, "GSEA_sym.txt", sep = "_"), row.names = FALSE, col.names = TRUE, sep = "\t")


###produce cls file which can be used in GSEA 

cls_1 <- paste(length(rt_M_S_l), 2, 1, sep = " ")
cls_2 <- paste(paste("#", cancer1, sep = ""), cancer2, sep = " ")
cls_3 <- paste(paste(rep(cancer1, length(rt_c1_T)), collapse = " "), paste(rep(cancer2, length(rt_c2_T)), collapse = " "), collapse = " ")
cls_file <- rbind(cls_1, cls_2, cls_3)
write.table(cls_file, file = paste(cancer1, cancer2, "group.cls", sep = "_"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  
  
  
  
  
  
  






















