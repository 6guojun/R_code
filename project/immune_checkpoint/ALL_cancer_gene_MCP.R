library(devtools)
library(MCPcounter)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/MCP")

cancer_list <- c("UVM", "GBM", "LUAD")

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
#base on the order of pdl1 pd-1 and CTLA4
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

####################################################
#correaltion plot betwwen socres and PDL1, PD1 and CTLA4
#
####################################################
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/Mut_correlation.R")
rt_score_genes <- rbind(rt_score_T, PDL1, PD1, CTLA4, FOXP3) 
MakCorPlot(t(rt_score_genes), Tnam = paste("cancer", "_gene_score_cor", sep = ""), height = 15, width = 15, cex.axis = 2)




###cluster heatmap of score
names(cancer_type) <- colnames(rt_score_T)
mat_score = as.matrix(rt_score_T)
mat_scaled = t(apply(mat_score, 1, scale))

ha1 = HeatmapAnnotation(cancer_type = cancer_type, PDL1 = PDL1 , PD1 = PD1 , FOXP3 = FOXP3 , CTLA4 = CTLA4 , 
                        col = list(
                          cancer_type = structure(names = unique(cancer_type ), colors()[2: (length(unique(cancer_type_PDL1)) + 1)]),
                          PDL1 = colorRamp2(c(0, max(PDL1_PDL1)/2), c("white", "#FFB400")),
                          PD1 = colorRamp2(c(0, max(PD1 )/2), c("white", "#007A87")),
                          FOXP3 = colorRamp2(c(0, max(FOXP3 )/2), c("white", "#FFAA91"))),
                        CTLA4 = colorRamp2(c(0, max(CTLA4 )/2), c("white", "#FF5A5F")))     


pdf(file = paste("cancer_cluster", "_KEGG_cluster.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   top_annotation = ha1, top_annotation_height = unit(2, "cm"), 
                   cluster_columns = TRUE, column_dend_reorder = TRUE, 
                   cluster_rows = TRUE, row_dend_reorder = TRUE, 
                   show_row_dend = TRUE, show_column_dend = TRUE,
                   show_row_names = TRUE, show_column_names = FALSE) 
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

####heatmap plot base on the expression of PDL1, PD1 and CTLA4
order_PDL1 <- order(PDL1)
order_PD1 <- order(PD1)
order_CTLA4 <- order(CTLA4)
order_FOXP3 <- order(FOXP3)

###PDL1 heatmap
rt_score_T_PDL1 <- rt_score_T[, order_PDL1]
cancer_type_PDL1 <- cancer_type[order_PDL1]
names(cancer_type_PDL1) <- colnames(rt_score_T_PDL1)
mat_PDL1 = as.matrix(rt_score_T_PDL1)
mat_scaled_PDL1 = apply(mat_PDL1, 1, scale)
PDL1_PDL1 <- PDL1[order_PDL1]
PD1_PDL1 <- PD1[order_PDL1]
CTLA4_PDL1 <- CTLA4[order_PDL1]
FOXP3_PDL1 <- FOXP3[order_PDL1]


mat_score_scale <- t(apply(rt_score_T_PDL1, 1, scale))

ha1 = HeatmapAnnotation(cancer_type = cancer_type_PDL1, PDL1 = PDL1_PDL1, PD1 = PD1_PDL1, FOXP3 = FOXP3_PDL1, CTLA4 = CTLA4_PDL1, 
                        col = list(
                          cancer_type = structure(names = unique(cancer_type_PDL1), colors()[2: (length(unique(cancer_type_PDL1)) + 1)]),
                          PDL1 = colorRamp2(c(0, max(PDL1_PDL1)/2), c("white", "#FFB400")),
                          PD1 = colorRamp2(c(0, max(PD1_PDL1)/2), c("white", "#007A87")),
                          FOXP3 = colorRamp2(c(0, max(FOXP3_PDL1)/2), c("white", "#FFAA91"))),
                        CTLA4 = colorRamp2(c(0, max(CTLA4_PDL1)/2), c("white", "#FF5A5F")))     


pdf(file = paste("cancer_PDL1", "_KEGG_cluster.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   top_annotation = ha1, top_annotation_height = unit(2, "cm"), 
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = FALSE,
                   show_row_names = TRUE, show_column_names = FALSE) 
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

###CTLA4 heatmap
rt_score_T_CTLA4 <- rt_score_T[, order_CTLA4]
cancer_type_CTLA4 <- cancer_type[order_CTLA4]
names(cancer_type_CTLA4) <- colnames(rt_score_T_CTLA4)
mat_CTLA4 = as.matrix(rt_score_T_CTLA4)
mat_scaled_CTLA4 = apply(mat_CTLA4, 1, scale)
PDL1_CTLA4 <- PDL1[order_CTLA4]
PD1_CTLA4 <- PD1[order_CTLA4]
CTLA4_CTLA4 <- CTLA4[order_CTLA4]
FOXP3_CTLA4 <- FOXP3[order_CTLA4]


mat_score_scale <- t(apply(rt_score_T_CTLA4, 1, scale))

ha1 = HeatmapAnnotation(cancer_type = cancer_type_CTLA4, PDL1 = PDL1_CTLA4, PD1 = PD1_CTLA4, FOXP3 = FOXP3_CTLA4, CTLA4 = CTLA4_CTLA4, 
                        col = list(
                          cancer_type = structure(names = unique(cancer_type_CTLA4), colors()[2: (length(unique(cancer_type_PDL1)) + 1)]),
                          PDL1 = colorRamp2(c(0, max(PDL1_PDL1)/2), c("white", "#FFB400")),
                          PD1 = colorRamp2(c(0, max(PD1_CTLA4)/2), c("white", "#007A87")),
                          FOXP3 = colorRamp2(c(0, max(FOXP3_CTLA4)/2), c("white", "#FFAA91"))),
                        CTLA4 = colorRamp2(c(0, max(CTLA4_CTLA4)/2), c("white", "#FF5A5F")))     


pdf(file = paste("cancer_CTLA4", "_KEGG_cluster.pdf", sep = ""), 10, 10)
ht_list <- Heatmap(mat_score_scale, name = "expression", col = colorRamp2(c(-5 , 0, 5), c(" blue", "white", "red")),
                   top_annotation = ha1, top_annotation_height = unit(2, "cm"), 
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = FALSE,
                   show_row_names = TRUE, show_column_names = FALSE) 
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()
