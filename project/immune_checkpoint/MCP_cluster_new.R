library(devtools)
library(MCPcounter)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/MCP")

cancer_list <- c("UVM", "GBM", "LUAD")
color_type_table <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/color.txt", sep = "\t", stringsAsFactors = FALSE)
color_type <- color_type_table$V2
####################################################
#get cancer exprssion 
#get the number of mutation and neo-ag
####################################################

####cancer expression and convert
source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")

GetBigMat <- function(cancer_list){
  cancer <- cancer_list[1]
  rt_on <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""), 
                      sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  rt_T <- split_tcga_tn(rt_on, sam_type = "tumor", split_type = "[.]")
  print(c(cancer, length(rt_T[1,])))
  rt_on_m <- data.frame(cbind(cancer, t(rt_T)), stringsAsFactors = FALSE)
  
  for(i in 2: length(cancer_list)){
    ###add the other cancer types
    cancer_add <- cancer_list[i]
    rt_add <-read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer_add, "/RNA_seq/TCGA-", cancer_add, ".htseq_fpkm.tsv", sep = ""), 
                        sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rt_add_T <- split_tcga_tn(rt_add, sam_type = "tumor", split_type = "[.]")
    print(c(cancer_add, length(rt_add_T[1, ])))
    rt_add_m <- data.frame(cbind(cancer_add, t(rt_add_T)), stringsAsFactors = FALSE)
    colnames(rt_add_m)[1] <- "cancer"
    if(all(colnames(rt_on_m) == colnames(rt_add_m))){
      rt_on_m <- data.frame(rbind(rt_on_m, rt_add_m), stringsAsFactors = FALSE)
    } else {
      stop('the gene names and order of all cancer must be same')
    }
  }
  return(rt_on_m)
}

rt_caner_all <- GetBigMat(cancer_list)

ensemble_gene <- matrix(unlist(strsplit(colnames(rt_caner_all)[-1], "[.]")), ncol = 2, byrow = TRUE)[, 1]
eg = bitr(ensemble_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
symbol_gene <- eg$SYMBOL[match(ensemble_gene, eg$ENSEMBL, nomatch = 0)]
colnames(rt_caner_all) <- c("group", paste(ensemble_gene, symbol_gene, sep = "_"))


###################################################
#get score matrix
#
###################################################
make_mcp <- function(rt_mcp){
  ###your gene ID should be set gene sumbol
  rt_score <- MCPcounter.estimate(rt_mcp,featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[2],
                                  probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                                  genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
  )
  return(rt_score)
}



###################################################
#MCP counter all samples included tumor 
#
###################################################
mcp_gene <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/MCP/MCP_counter_signature.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mcp_gene_pos <- match(mcp_gene$HUGO.symbols , symbol_gene, nomatch = 0)
                      
###the rownames of rt_T/N is the same with rt_sym
rt_T_cmp <- rt_caner_all[, -1][, mcp_gene_pos]
colnames(rt_T_cmp) <- symbol_gene[mcp_gene_pos]
#eg_test <- bitr(matrix(unlist(strsplit(colnames(rt_T_cmp), "[.]")), ncol = 2, byrow = TRUE)[, 1], fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
#rt_N_cmp <- rt_N[mcp_gene_pos, ]
#row.names(rt_N_cmp) <- rt_sym$SYMBOL[mcp_gene_pos]

#rt_score_N <- make_mcp(rt_N_cmp)

cancer_type <- rt_caner_all$group  ###colnames of rt_score_T share the same samples names with rownames of rt_cancer_all
mat_mcp <- apply(rt_T_cmp, 1, function(x){as.numeric(x)})
rt_T_cmp_V <- data.frame(mat_mcp, stringsAsFactors = FALSE)
row.names(rt_T_cmp_V) <- colnames(rt_T_cmp)
rt_score_T <- make_mcp(rt_T_cmp_V)

#tumor
pdf(file = paste("cancer", "cmp_sore_tumor_raw_2.pdf", sep = "_"), 10, 5)
heatmap(as.matrix(rt_score_T),col=colorRampPalette(c("blue","white","red"))(100))
dev.off()


###################################################
#get expression of pdl1 pd-1 and CTLA4
#PD1 is ENSG00000188389_PDCD1
#
###################################################
PDL1 <- as.numeric(rt_caner_all[, grep("CD274", colnames(rt_caner_all))])
PD1  <- as.numeric(rt_caner_all[, grep("ENSG00000188389", colnames(rt_caner_all))])
CTLA4 <- as.numeric(rt_caner_all[, grep("CTLA4", colnames(rt_caner_all))])
FOXP3 <- as.numeric(rt_caner_all[, grep("FOXP3", colnames(rt_caner_all))])


names(PDL1) <- row.names(rt_caner_all)
names(PD1) <- row.names(rt_caner_all)
names(CTLA4) <- row.names(rt_caner_all)
names(FOXP3) <- row.names(rt_caner_all)

##################################################
#MCP score mean or median of all cancer
#
##################################################
#colnames of rt_score_T is the same with cancer type , so we can get each cancer sore from rt_score based on cancer type
CountCmeanScore <- function(cancer_list, cancer_type, rt_score_T){
  cancer_one <- cancer_list[1]
  print(cancer_one)
  cancer_pos_one <- grep(cancer_one, cancer_type)
  cancer_score_mean_one <- apply(rt_score_T[, cancer_pos_one], 1, mean)
  for(i in 2: length(cancer_list)){
    ###add the other cancer types
    cancer_oth <- cancer_list[i]
    print(cancer_oth)
    cancer_pos_oth <- grep(cancer_oth, cancer_type)
    cancer_score_mean_oth <- apply(rt_score_T[, cancer_pos_oth], 1, mean)
    cancer_score_mean_one <- data.frame(cbind(cancer_score_mean_one, cancer_score_mean_oth), stringsAsFactors = FALSE)
  }
  colnames(cancer_score_mean_one) <- cancer_list
  return(cancer_score_mean_one)
}
  
cancer_score_mean <- CountCmeanScore(cancer_list, cancer_type, rt_score_T)
mat_score_scale_all <- t(apply(cancer_score_mean, 1, scale))

###the name of PDL1,CTLA4 is the same name with the rowname of rt_cancer_all
VecTolist <- function(cancer_list, cancer_type, gene_exp, out_type = c("exp", "median")){
   one_cancer_pos <- grep(cancer_list[1], cancer_type)
   exp_cancer_one <- gene_exp[one_cancer_pos]
   exp_list <- list(exp_cancer_one)
   names(exp_list) <- cancer_list[1]
   #count median
   exp_median <- median(exp_cancer_one)
   names(exp_median) <- cancer_list[1]
  for(i in 2: length(cancer_list)){
     oth_cancer_pos <- grep(cancer_list[i], cancer_type)
     exp_cancer_oth <- gene_exp[oth_cancer_pos]
     exp_list_oth <- list(exp_cancer_oth)
     names(exp_list_oth) <- cancer_list[i]
     exp_list <- c(exp_list, exp_list_oth)
     #count median
     exp_median_oth <- median(exp_cancer_oth)
     names(exp_median_oth) <- cancer_list[i]
     exp_median <- c(exp_median, exp_median_oth)
  }
   if(out_type == "exp"){
     return(exp_list)
   } else {
     return(exp_median)
   }
}

CTLA4_list <- VecTolist(cancer_list, cancer_type, CTLA4, out_type = "exp")
CTLA4_order <- order(VecTolist(cancer_list, cancer_type, CTLA4, out_type= "median"))
PDL1_list <- VecTolist(cancer_list, cancer_type, PDL1)
PDL1_order <- order(VecTolist(cancer_list, cancer_type, PDL1, out_type= "median"))
PD1_list <- VecTolist(cancer_list, cancer_type, PD1)
PD1_order <- order(VecTolist(cancer_list, cancer_type, PD1, out_type= "median"))


ha1 = HeatmapAnnotation(CTLA4 = anno_boxplot(CTLA4_list, axis = TRUE, gp = gpar(fill= color_type[1: (length(unique(cancer_list)))])), 
                        PDL1 = anno_boxplot(PDL1_list, axis = TRUE, gp = gpar(fill= color_type[1: (length(unique(cancer_list)))])), 
                        PD1 = anno_boxplot(PD1_list, axis = TRUE, gp = gpar(fill = color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)

ha2 = HeatmapAnnotation(cancer_name = cancer_list, 
                        col = list(cancer_name = structure(names = cancer_list, color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)                           

pdf(file = paste("cancer_mean", "_MCPscore_cluster.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale_all, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   bottom_annotation = ha1, bottom_annotation_height = unit(10, "cm"), 
                   top_annotation = ha2, top_annotation_height = unit(0.5, "cm"),
                   cluster_columns = TRUE, column_dend_reorder = TRUE, 
                   cluster_rows = TRUE, row_dend_reorder = TRUE, 
                   show_row_dend = TRUE, show_column_dend = TRUE,
                   show_row_names = TRUE, show_column_names = TRUE)
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

###base on PDL1 order 
PDL1_list_PDL1 <- PDL1_list[PDL1_order]
CTLA4_list_PDL1 <- CTLA4_list[PDL1_order]
PD1_list_PDL1 <- PD1_list[PDL1_order]
cancer_score_mean_DPL1 <- CountCmeanScore(cancer_list, cancer_type, rt_score_T)[PDL1_order]
mat_score_scale_all_PDL1 <- t(apply(cancer_score_mean_DPL1, 1, scale))

###make plot
ha1 = HeatmapAnnotation(PDL1 = anno_boxplot(PDL1_list_PDL1, axis = TRUE, gp = gpar(col= color_type[1: (length(unique(cancer_list)))])), 
                        CTLA4 = anno_boxplot(CTLA4_list_PDL1, axis = TRUE, gp = gpar(col= color_type[1: (length(unique(cancer_list)))])), 
                        PD1 = anno_boxplot(PD1_list_PDL1, axis = TRUE, gp = gpar(col = color_type[1: (length(unique(cancer_list)))])), 
                        cancer_type = anno_text(names(PDL1_list_PDL1), rot = 90, just = "left", offset = unit(2, "mm")), 
                        show_annotation_name = TRUE)
cancer_list_PDL1 <- cancer_list[PDL1_order]
ha2 = HeatmapAnnotation(cancer_name = cancer_list_PDL1, 
                        col = list(cancer_name = structure(names = cancer_list_PDL1, color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)                           

pdf(file = paste("cancer_mean", "_MCPscore_PDL1.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale_all_PDL1, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   bottom_annotation = ha1, bottom_annotation_height = unit(12, "cm"), 
                   top_annotation = ha2, top_annotation_height = unit(0.5, "cm"),
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = TRUE,
                   show_row_names = TRUE, show_column_names = TRUE)
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

###base on PD1 order 
PDL1_list_PD1 <- PDL1_list[PD1_order]
CTLA4_list_PD1 <- CTLA4_list[PD1_order]
PD1_list_PD1 <- PD1_list[PD1_order]
cancer_score_mean_DPL1 <- CountCmeanScore(cancer_list, cancer_type, rt_score_T)[PD1_order]
mat_score_scale_all_PD1 <- t(apply(cancer_score_mean_DPL1, 1, scale))

###make plot
ha1 = HeatmapAnnotation(PD1 = anno_boxplot(PD1_list_PD1, axis = TRUE, gp = gpar(col = color_type[1: (length(unique(cancer_list)))])), 
                        PDL1 = anno_boxplot(PDL1_list_PD1, axis = TRUE, gp = gpar(col= color_type[1: (length(unique(cancer_list)))])), 
                        CTLA4 = anno_boxplot(CTLA4_list_PD1, axis = TRUE, gp = gpar(col= color_type[1: (length(unique(cancer_list)))])), 
                        cancer_type = anno_text(names(PD1_list_PD1), rot = 90, just = "left", offset = unit(2, "mm")), 
                        show_annotation_name = TRUE)

cancer_list_PD1 <- cancer_list[PD1_order]
ha2 = HeatmapAnnotation(cancer_name = cancer_list_PD1, 
                        col = list(cancer_name = structure(names = cancer_list_PD1, color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)                           

pdf(file = paste("cancer_mean", "_MCPscore_PD1.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale_all_PD1, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   bottom_annotation = ha1, bottom_annotation_height = unit(12, "cm"), 
                   top_annotation = ha2, top_annotation_height = unit(0.5, "cm"),
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = TRUE,
                   show_row_names = TRUE, show_column_names = TRUE)
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

###base on CTLA4 order 
PDL1_list_CTLA4 <- PDL1_list[CTLA4_order]
CTLA4_list_CTLA4 <- CTLA4_list[CTLA4_order]
PD1_list_CTLA4 <- PD1_list[CTLA4_order]
cancer_score_mean_DPL1 <- CountCmeanScore(cancer_list, cancer_type, rt_score_T)[CTLA4_order]
mat_score_scale_all_CTLA4 <- t(apply(cancer_score_mean_DPL1, 1, scale))

###make plot
ha1 = HeatmapAnnotation(CTLA4 = anno_boxplot(CTLA4_list_CTLA4, axis = TRUE, gp = gpar(col= color_type[1: (length(unique(cancer_list)))])), 
                        PD1 = anno_boxplot(PD1_list_CTLA4, axis = TRUE, gp = gpar(col = color_type[1: (length(unique(cancer_list)))])), 
                        PDL1 = anno_boxplot(PDL1_list_CTLA4, axis = TRUE, gp = gpar(col= color_type[1: (length(unique(cancer_list)))])), 
                        cancer_type = anno_text(names(CTLA4_list_CTLA4), rot = 90, just = "left", offset = unit(2, "mm")), 
                        show_annotation_name = TRUE)

cancer_list_CTLA4 <- cancer_list[CTLA4_order]
ha2 = HeatmapAnnotation(cancer_name = cancer_list_CTLA4, 
                        col = list(cancer_name = structure(names = cancer_list_CTLA4, color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)      

pdf(file = paste("cancer_mean", "_MCPscore_CTLA4.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale_all_CTLA4, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   bottom_annotation = ha1, bottom_annotation_height = unit(12, "cm"), 
                   top_annotation = ha2, top_annotation_height = unit(0.5, "cm"),
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = TRUE,
                   show_row_names = TRUE, show_column_names = TRUE)
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()


