setwd("/Users/stead/Desktop/GSEA")
cancer <- "GBM"
rt <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                 sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt$Ensembl_ID <- matrix(unlist(strsplit(rt$Ensembl_ID, "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(rt$Ensembl_ID, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
eg_rm <- eg[!is.na(eg$SYMBOL), ]

###get intersection genes of gene symbols
com_genes <- intersect(eg_rm$ENSEMBL, rt$Ensembl_ID)
eg_com <- eg_rm[match(com_genes, eg_rm$ENSEMBL, nomatch = 0), ]
rt_com <- rt[match(com_genes, rt$Ensembl_ID, nomatch = 0), -1]


###reorder tumor and normal samples to perform GSEA analysis
get_sam_nam <- function(rt_tcga, split_type){
  if(split_type == "[.]"){
    sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "[.]")), ncol = 4, byrow = TRUE)
  }else if(split_type == "-"){
    sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "-")), ncol = 4, byrow = TRUE)
  } else if(split_type == "_") {
    sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "_")), ncol = 4, byrow = TRUE)
  } else {
    stop("you must set split_type ")
  }
  sam_nam <- paste(sam_G[, 1], sam_G[, 2], sam_G[, 3], sep = "-")
  sam_list <- list(sam_G, sam_nam)
  return(sam_list)
}


sam_list <- get_sam_nam(rt_com, "[.]")
sam_G <- sam_list[[1]]

tcga_pos_T <- c(grep("01", sam_G[, 4]), grep("02", sam_G[, 4]), grep("03", sam_G[, 4]), grep("04", sam_G[, 4]), grep("05", sam_G[, 4]), 
                grep("06", sam_G[, 4]), grep("07", sam_G[, 4]), grep("08", sam_G[, 4]), grep("09", sam_G[, 4]))
tcga_pos_N <- c(grep("11", sam_G[, 4]), grep("10", sam_G[, 4]))

rt_GSEA <- rt_com[, c(tcga_pos_N, tcga_pos_T)]
rt_sym_GSEA <- cbind(eg_com$SYMBOL, rt_GSEA)
colnames(rt_sym_GSEA)[1] <- "gene_id"
write.table(rt_sym_GSEA, file = paste(cancer, "GSEA_sym.txt", sep = "_"), row.names = FALSE, col.names = TRUE, sep = "\t")

###produce cls file which can be used in GSEA 
cls_1 <- paste(length(rt_sym_GSEA[, -1]), 2, 1, sep = " ")
cls_2 <- paste("#normal", "tumor", sep = " ")
cls_3 <- paste(paste(rep("normal", length(tcga_pos_N)), collapse = " "), paste(rep("tumor", length(tcga_pos_T)), collapse = " "), collapse = " ")
cls_file <- rbind(cls_1, cls_2, cls_3)
write.table(cls_file, file = paste(cancer, "_group.cls", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

