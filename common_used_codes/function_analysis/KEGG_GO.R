###you shold prepare a fc and pvalue table which contain gene name, pval and FC 
###you must convert your geneID to ENTREZID before you do kegg and go analysis
###packages:clusterProfiler, org.Hs.eg.db, mygene, annotables, pathview
###ENTREZID: 4312
###ENSEMBLEID: ENSG00000211445
###ENSEMBL: GPX3
###Usemethod
###MakKEGO(rt_KG, ont = ont, data_type = "gene_nam", Gene_ID = "ENSEMBL", OrgDb = "org.Rn.eg.db)
###MakKEGO(rt_KG, ont = ont, data_type = "gene_nam", Gene_ID = "ENSEMBL", OrgDb = "org.Rn.eg.db)
###MakKEGO(rt_KG, ont = ont, data_type = "gene_list", Gene_ID = "SYMBOL", OrgDb = "org.Rn.eg.db)
###species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
###reference: https://bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

library(clusterProfiler)
#library(org.Rn.eg.db)
library(mygene)
library(annotables)
library(pathview)

source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_coord.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_E.R")
########################################################################
#rt <- read.table(file = "/Users/stead/Desktop/LungCancer_SS/LUAD/DEG/WT_Mut_DEG_Pval_FC.txt", 
#                 header = TRUE, row.names = 1, sep = "\t")
#row.names(rt) <- matrix(unlist(strsplit(row.names(rt), "\\.")), ncol = 2, byrow = TRUE)[, 1]
#######################################################################

##########################KEGG and GO analysis###########################
MakKEGO <- function(data_type = c("gene_nam", "gene_list"), gene_fun = gene_fun, genes_list_fun = genes_list_fun, num_term = 10, main_tile = NULL){
  #if you don't have gene_list which contain FC or log2FC (we need the colname is log2FC), you can set rt_KG as gene names
  #so you will get the results about gene names
  #if you have gene_list wihch contain log2FC, you can get all results
  #ont is one of 'MF', 'BP', 'CC' subontologies 
  
  kk <- enrichKEGG(gene = gene_fun, organism = 'hsa', pvalueCutoff = 0.05)
  
  mat_kk_all <- data.frame(cbind(kk@result$Description, kk@result$Count, kk@result$p.adjust, kk@result$pvalue), stringsAsFactors = FALSE)
  colnames(mat_kk_all) <- c('KEGG_term', 'genes', "p_adjust", 'pvalue')
  write.table(mat_kk_all, file = 'kk_res.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
  mat_kk <- mat_kk_all[1:num_term, ]
  colnames(mat_kk) <- c('KEGG_term', 'genes', "p_adjust", 'pvalue')
  mat_kk$genes <- as.numeric(mat_kk$genes)
  mat_kk$p_adjust <- as.numeric(format(as.numeric(mat_kk$p_adjust), digits = 2))
  mat_kk$pvalue <- as.numeric(format(as.numeric(mat_kk$pvalue), digits = 2))
  
  ggo_CC <- groupGO(gene = gene_fun, OrgDb = org.Hs.eg.db, ont = "CC", level = 3, readable = TRUE)
  ggo_BP <- groupGO(gene = gene_fun, OrgDb = org.Hs.eg.db, ont = "BP", level = 3, readable = TRUE)
  ggo_MF <- groupGO(gene = gene_fun, OrgDb = org.Hs.eg.db, ont = "MF", level = 3, readable = TRUE)
  
  ggo_CBM <-  data.frame(rbind(cbind('CC', as.character(ggo_CC@result$Description), ggo_CC@result$Count), 
                               cbind('BP', as.character(ggo_BP@result$Description), ggo_BP@result$Count), 
                               cbind('MF', as.character(ggo_MF@result$Description), ggo_MF@result$Count)), stringsAsFactors = FALSE)
  colnames(ggo_CBM) <- c('GO_type', 'GO_term', "genes")
  write.table(ggo_CBM, file = 'CC_BP_MF_ggo.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
  
  ggo_all <-  data.frame(rbind(cbind('CC', as.character(ggo_CC@result$Description), ggo_CC@result$Count)[1:num_term, ], 
                               cbind('BP', as.character(ggo_BP@result$Description), ggo_BP@result$Count)[1:num_term, ], 
                               cbind('MF', as.character(ggo_MF@result$Description), ggo_MF@result$Count)[1:num_term, ]), stringsAsFactors = FALSE)
  
  colnames(ggo_all) <- c('GO_type', 'GO_term', "genes")
  ggo_all$genes <- as.numeric(ggo_all$genes)
  
  pdf(file = 'KEGG_bar_pvalue.pdf', 10, 5)
  kegg_bar_pv <- ggplot(mat_kk, aes(x = reorder(KEGG_term, pvalue), y = genes, fill = pvalue)) + geom_bar(stat = "identity") + 
    xlab('KEGG_term') + ylab('genes') + ggtitle(main_tile) + 
    coord_flip() + theme_coord 
  print(kegg_bar_pv)
  dev.off()
  
  pdf(file = 'KEGG_bar_padjust.pdf', 10, 5)
  kegg_bar_pj <- ggplot(mat_kk, aes(x = reorder(KEGG_term, p_adjust), y = genes, fill = p_adjust)) + geom_bar(stat = "identity") + 
    xlab('KEGG_term') + ylab('genes') +  ggtitle(main_tile) + 
    coord_flip() + theme_coord 
  print(kegg_bar_pj)
  dev.off()
  
  pdf(file = 'ggo_bar.pdf', 10, 5)
  ggo_bar <-  ggplot(ggo_all, aes(x = GO_term, y = genes, fill = GO_type)) + geom_bar(stat = "identity") + 
    facet_wrap(~ GO_type, scales = "free") +  xlab('GO_term') + ylab('genes') +  ggtitle(main_tile) + 
    theme_classic() + theme(axis.title.x = element_text(size = 5, color = "black", vjust = 0.5, hjust = 0.5, angle = 1),
                            axis.text.x = element_text(size = 5, color = "black", vjust = 0.5, hjust = 0.5, angle = 90)) 
   print(ggo_bar)
  dev.off()
  
  if(data_type == "gene_nam"){
    print('results finish')
    
  } 
  else if (data_type == "gene_list"){
    
    kk <- enrichKEGG(gene = gene_fun, organism = 'hsa', pvalueCutoff = 0.05)
    ggo <- groupGO(gene = gene_fun, OrgDb = org.Hs.eg.db, ont = "CC", level = 3, readable = TRUE)
    
    ego_CC <- enrichGO(gene = gene_fun, universe = names(genes_list_fun), OrgDb = org.Hs.eg.db, ont = "CC",
                       pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1, readable = TRUE)
    ego_BP <- enrichGO(gene = gene_fun, universe = names(genes_list_fun), OrgDb = org.Hs.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1, readable = TRUE)
    ego_MF <- enrichGO(gene = gene_fun, universe = names(genes_list_fun), OrgDb = org.Hs.eg.db, ont = "MF",
                       pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.1, readable = TRUE)
    
    ego_CBM <-  data.frame(rbind(cbind('CC', ego_CC@result$Description, ego_CC@result$Count, ego_CC@result$p.adjust, ego_CC@result$pvalue), 
                                 cbind('BP', ego_BP@result$Description, ego_BP@result$Count, ego_BP@result$p.adjust, ego_BP@result$pvalue), 
                                 cbind('MF', ego_MF@result$Description, ego_MF@result$Count, ego_MF@result$p.adjust, ego_MF@result$pvalue)), stringsAsFactors = FALSE)
    colnames(ego_CBM) <- c('GO_type', 'GO_term', "Genes", "p_adjust", "p_value") 
    write.table(ego_CBM, file = 'CC_BP_MF_ego.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
    
    ego_all <-  data.frame(rbind(cbind('CC', ego_CC@result$Description, ego_CC@result$Count, ego_CC@result$p.adjust, ego_CC@result$pvalue)[1:num_term, ], 
                                 cbind('BP', ego_BP@result$Description, ego_BP@result$Count, ego_BP@result$p.adjust, ego_BP@result$pvalue)[1:num_term, ], 
                                 cbind('MF', ego_MF@result$Description, ego_MF@result$Count, ego_MF@result$p.adjust, ego_MF@result$pvalue)[1:num_term, ]), stringsAsFactors = FALSE)
    colnames(ego_all) <- c('GO_type', 'GO_term', "Genes", "p_adjust", "p_value") 
    ego_all$Genes <- as.numeric(ego_all$Genes)
    ego_all$p_adjust <- as.numeric(format(as.numeric(ego_all$p_adjust), digits = 2))
    ego_all$p_value <- as.numeric(format(as.numeric(ego_all$p_value), digits = 2))
    
    pdf(file = 'ego_bar.pdf')
    ego_bar <-  ggplot(ego_all, aes(x = GO_term, y = Genes, fill = p_adjust)) + geom_bar(stat = "identity") +
      facet_wrap(~ GO_type, scales = "free") + 
      theme_classic() + theme(axis.title.x = element_text(size = 5, color = "black", vjust = 0.5, hjust = 0.5, angle = 1),
                              axis.text.x = element_text(size = 5, color = "black", vjust = 0.5, hjust = 0.5, angle = 90))
    print(ego_bar)
    dev.off()
    
    pdf(file = 'CC_emapplot.pdf')
    em_ego_cc <-  emapplot(ego_CC)
    print(em_ego_cc)
    dev.off()
    
    pdf(file = 'BP_emapplot.pdf')
    em_ego_bp <- emapplot(ego_BP)
    print(em_ego_bp)
    dev.off()
    
    pdf(file = 'MF_emapplot.pdf')
    em_ego_mf <- emapplot(ego_MF)
    print(em_ego_mf)
    dev.off()
    
    pdf(file = 'CC_cnetplot.pdf')
    cen_cc <-  cnetplot(ego_CC, categorySize= "pvalue", foldChange = genes_list_fun)
    print(cen_cc)
    dev.off()
    
    pdf(file = 'BP_cnetplot.pdf')
    cen_bp <- cnetplot(ego_BP, categorySize= "pvalue", foldChange = genes_list_fun)
    print(cen_bp)
    dev.off()
    
    pdf(file = 'MF_cnetplot.pdf')
    cen_mf <- cnetplot(ego_MF, categorySize= "pvalue", foldChange = genes_list_fun)
    print(cen_mf)
    dev.off()
    
  }else {
    print('please set data_type')
  }
  }

####################################################################

