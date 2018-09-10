library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(annotables)
library(pathview)
library(RColorBrewer)
library(scales)

setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/correlation_and_go_cluster")

cancer_list <- c("GBM", "UVM") 
color_type_table <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/color.txt", sep = "\t", stringsAsFactors = FALSE)
color_type <- color_type_table$V2
immune_go <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/immune_pathway/go_term.txt", sep = "\t", stringsAsFactors = FALSE)
####################################################
#get a dataframe which comtain all tumor samples
#**the first column is ENTREZID**, colnames is samples id
###################################################
###usage: rt_one <- GetOneData(cancer)
source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_tcga_sample_type.R")
GetOneData <- function(cancer){
  rt<- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""), 
                  sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  rt_T <- split_tcga_tn(rt, sam_type = "tumor", split_type = "[.]")
  row.names(rt_T) <- matrix(unlist(strsplit(row.names(rt_T), "[.]")), ncol = 2, byrow = TRUE)[, 1]
  eg = bitr(row.names(rt_T), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  rt_ent <- merge(eg, cbind(row.names(rt_T), rt_T), by = 1)[, -1] 
  print(c(cancer, dim(rt_ent)))
  return(rt_ent)
}

#rt_one <- GetOneData("UVM")

####################################################
#obtain high correaltion gene and correlation value from a expression table
#
###################################################
###usage:get_high_cor_gene(in_gene, gene_exp, in_mat, cutoff)
get_high_cor_gene <- function(in_gene, gene_exp, in_mat, cutoff){
  #in_gene is an gene name 
  #gene_exp is one gene expression value such as CTLA4 expression
  #im_mat is an expression table which contian in_gene or other genes
  #**the first column of in_mat is SYMBOL**
  #cutoff was used to obtain genes which meet the request
  gene_exp <- as.numeric(gene_exp)
  test_gene_exp <- as.numeric(in_mat[in_gene, -1])
  if(all(test_gene_exp == 0)){
    return(NULL)
  }
  if(any(test_gene_exp != 0)){
    cor_value = round(cor(gene_exp, test_gene_exp), 2)
    if(cor_value > cutoff){
      gene_value <- c(in_gene, cor_value)
      return (gene_value)
    } else {
      gene_value <- c(in_gene, cor_value)
      return(NULL)
    } 
  }
}

#gene_value_one <- get_high_cor_gene(in_gene, gene_exp, in_mat, cutoff)
####################################################
#obtain immune related pathway 
#
###################################################
###usage:GetEnrich(cancer, in_mat, gene_exp, immune_go, cutoff = cutoff)
GetEnrich <- function(cancer, in_mat, gene_exp, immune_go, cutoff = cutoff){
  #cancer is one cancer name
  #rt_go contain pathway name and gene coutns which enrich in immune realted pathway
  in_gene = row.names(in_mat)
  cor_gene_list <- lapply(in_gene, get_high_cor_gene, gene_exp, in_mat, cutoff)
  rt_cor_gene <- do.call(rbind, cor_gene_list)
  gene <- in_mat[, 1][as.numeric(rt_cor_gene[, 1])]
  gene_ent <- bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  go <- groupGO(gene = gene_ent$ENTREZID, OrgDb = org.Hs.eg.db, keytype = "ENTREZID", on = "BP", level = 3, readable = TRUE)
  print(head(go))
  path_M <- intersect(go@result$ID, immune_go$V1)
  rt_go <- go@result[, c("ID", "Count")][match(path_M, go@result$ID, nomatch = 0), ] 
  colnames(rt_go) <- c("ID", cancer)
  return(rt_go)
}

#GetEnrich(cancer, in_mat, gene_exp, immune_go, cutoff = cutoff)

##################################################################
#main function which was used to get gene expression of all cancer 
#and get all immnue gene counts of all cancer
##################################################################
###usage:MainMergeBig(cancer_list, gene_SYMBOL, immune_go)
MainMergeBig <- function(cancer_list, gene_SYMBOL, immune_go){
  #cancer_list are many cancer names
  #gene_SYMBOL is one gene SYMBOL
  #immume_path is immune related pathway table
  cancer <- cancer_list[1]
  rt_ent_on <- GetOneData(cancer)
  gene_exp_on <- as.numeric(rt_ent_on[which(rt_ent_on$SYMBOL == gene_SYMBOL), -1])
  names(gene_exp_on) <- colnames(rt_ent_on[, -1])
  gene_exp_list <- list(gene_exp_on)
  names(gene_exp_list) <- paste(cancer, gene_SYMBOL, sep = "_")
  
  ##get immune related pathway gene count 
  rt_go_on <- GetEnrich(cancer, rt_ent_on, gene_exp_on, immune_go, cutoff = 0.4)
  rt_go_on_immune <- data.frame(merge(immune_go, rt_go_on, by = 1, all = TRUE), stringsAsFactors = FALSE)
  for(i in 2: length(cancer_list)){
    ###add the other cancer types
    cancer_add <- cancer_list[i]
    rt_ent_add <- GetOneData(cancer_add)
    gene_exp_add <- as.numeric(rt_ent_add[which( rt_ent_add$SYMBOL == gene_SYMBOL), -1])
    names(gene_exp_add) <- colnames(rt_ent_add[, -1])
    gene_exp_add_list <- list(gene_exp_add)
    names(gene_exp_add_list) <- paste(cancer_add, gene_SYMBOL, sep = "_")
    gene_exp_list <- c(gene_exp_list, gene_exp_add_list)
    
    ##get immune related pathway gene count 
    rt_go_add <- GetEnrich(cancer_add, rt_ent_add, gene_exp_add, immune_go, cutoff = 0.4)
    rt_go_on_immune <- data.frame(merge(rt_go_on_immune, rt_go_add, by = 1, all = TRUE), stringsAsFactors = FALSE)
  }
  immune_gene_mat_list <- list(gene_exp_list, rt_go_on_immune)
  return(immune_gene_mat_list)
}




#############################################################
#get gene expression 
#CTLA4 entrezid is 1493
############################################################
immune_gene_mat_list_CTLA4 <- MainMergeBig(cancer_list, "CTLA4", immune_go)
CTLA4_cancer_list <- immune_gene_mat_list_CTLA4[[1]]


#############################################################
#count immune realted pathway
#pathway results
############################################################
rt_immune_count_CTLA4 <- do.call(cbind, immune_gene_mat_list_CTLA4[2])
mat_immnue_CTLA4 <- rt_immune_count_CTLA4[, -c(1, 2)]
mat_immnue_CTLA4 <- apply(mat_immnue_CTLA4, 2, function(x){as.numeric(x)})
mat_immnue_CTLA4[is.na(mat_immnue_CTLA4)] <- 0
row.names(mat_immnue_CTLA4) <- rt_immune_count_CTLA4$V2

ha1 = HeatmapAnnotation(dist1 = anno_barplot(colSums(mat_immnue_CTLA4), bar_width = 1, gp = gpar(col = NA, fill = "skyblue"), 
                                             border = FALSE, axis = TRUE))

ha2 = HeatmapAnnotation(CTLA4 = anno_boxplot(CTLA4_cancer_list, axis = TRUE, gp = gpar(fill= color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)

pdf(file = paste("CTLA4", "_KEGG_count_cluster.pdf", sep = ""), 15, 10)
ht_list_CTLA4 <- Heatmap(mat_immnue_CTLA4, name = "gene count", col = colorRamp2(c(0, 20, 40, 60), c("grey80", "lightskyblue1", "steelblue", "firebrick2")),
                         cluster_columns = TRUE, show_row_dend = TRUE, rect_gp = gpar(col= "white"), show_column_names = TRUE, show_column_dend = FALSE, 
                         row_names_side = "right", row_names_gp = gpar(fontsize = 15), row_names_max_width = unit(18, "cm"),
                         top_annotation = ha1, top_annotation_height = unit(4, "cm"),
                         bottom_annotation = ha2, bottom_annotation_height = unit(5, "cm"),
                         column_title = 'number of immune pathway genes')

draw(ht_list_CTLA4, annotation_legend_side = "right", heatmap_legend_side = "left")
dev.off()


################################################################
#perform PDCD1 ('5133')
#
################################################################
immune_gene_mat_list_PD1 <- MainMergeBig(cancer_list, "PDCD1", immune_go)
PD1_cancer_list <- immune_gene_mat_list_PD1[[1]]

#############################################################
#count immune realted pathway
#pathway results
############################################################
rt_immune_count_PD1 <- do.call(cbind, immune_gene_mat_list_PD1[2])
mat_immnue_PD1 <- rt_immune_count_PD1[, -c(1, 2)]
mat_immnue_PD1 <- apply(mat_immnue_PD1, 2, function(x){as.numeric(x)})
mat_immnue_PD1[is.na(mat_immnue_PD1)] <- 0
row.names(mat_immnue_PD1) <- rt_immune_count_PD1$V2

ha1 = HeatmapAnnotation(dist1 = anno_barplot(colSums(mat_immnue_PD1), bar_width = 1, gp = gpar(col = NA, fill = "lightblue"), 
                                             border = FALSE, axis = TRUE))

ha2 = HeatmapAnnotation(PD1 = anno_boxplot(PD1_cancer_list, axis = TRUE, gp = gpar(fill= color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)

pdf(file = paste("PD1", "_KEGG_count_cluster.pdf", sep = ""), 15, 10)
ht_list_PD1 <- Heatmap(mat_immnue_PD1, name = "gene count", col = colorRamp2(c(0, 20, 40, 60), c("grey80", "lightskyblue1", "steelblue", "firebrick2")),
                       cluster_columns = TRUE, show_row_dend = TRUE, rect_gp = gpar(col= "white"), show_column_names = TRUE, show_column_dend = FALSE, 
                       row_names_side = "right", row_names_gp = gpar(fontsize = 15), row_names_max_width = unit(18, "cm"),
                       top_annotation = ha1, top_annotation_height = unit(1.5, "cm"),
                       bottom_annotation = ha2, bottom_annotation_height = unit(3, "cm"),
                       column_title = 'number of immune pathway genes')

draw(ht_list_PD1, annotation_legend_side = "right", heatmap_legend_side = "left")
dev.off()


################################################################
#perform CD274 ('29126')
#
################################################################
immune_gene_mat_list_PDL1 <- MainMergeBig(cancer_list, "CD274", immune_go)
PDL1_cancer_list <- immune_gene_mat_list_PDL1[[1]]


#############################################################
#count immune realted pathway
#pathway results
############################################################
rt_immune_count_PDL1 <- do.call(cbind, immune_gene_mat_list_PDL1[2])
mat_immnue_PDL1 <- rt_immune_count_PDL1[, -c(1, 2)]
mat_immnue_PDL1 <- apply(mat_immnue_PDL1, 2, function(x){as.numeric(x)})
mat_immnue_PDL1[is.na(mat_immnue_PDL1)] <- 0
row.names(mat_immnue_PDL1) <- rt_immune_count_PDL1$V2

ha1 = HeatmapAnnotation(dist1 = anno_barplot(colSums(mat_immnue_PDL1), bar_width = 1, gp = gpar(col = NA, fill = "lightblue"), 
                                             border = FALSE, axis = TRUE))

ha2 = HeatmapAnnotation(PDL1 = anno_boxplot(PDL1_cancer_list, axis = TRUE, gp = gpar(fill= color_type[1: (length(unique(cancer_list)))])), 
                        show_annotation_name = TRUE)

pdf(file = paste("PDL1", "_KEGG_count_cluster.pdf", sep = ""), 15, 10)
ht_list_PDL1 <- Heatmap(mat_immnue_PDL1, name = "gene count", col = colorRamp2(c(0, 10, 40, 80), c("grey80", "cornflowerblue", "skyblue", "red")),
                        cluster_columns = TRUE, show_row_dend = TRUE, rect_gp = gpar(col= "white"), show_column_names = TRUE, show_column_dend = FALSE, 
                        row_names_side = "right", row_names_gp = gpar(fontsize = 15), row_names_max_width = unit(18, "cm"),
                        top_annotation = ha1, top_annotation_height = unit(1.5, "cm"),
                        bottom_annotation = ha2, bottom_annotation_height = unit(3, "cm"),
                        column_title = 'number of immune pathway genes')

draw(ht_list_PDL1, annotation_legend_side = "right", heatmap_legend_side = "left")
dev.off()