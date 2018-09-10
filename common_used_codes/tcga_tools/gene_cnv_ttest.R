###do mut complex heatmap contained mutation, clinial and cnv
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_B.R")
source('/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/count_gene_mutation_rate.R')

PloCNVRisk <- function(cacner, rt_risk_cnv_order){
  rt_cnv_gene_c <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/', cancer, '/copy_num/Gistic2_CopyNumber_Gistic2_all_data_by_genes', sep = ""),  
                              header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
  rt_cnv_gene_l <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/',cancer, '/copy_num/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes', sep = ''), 
                              header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
  
  cnv_gene_c_id <- matrix(unlist(strsplit(colnames(rt_cnv_gene_c), '[.]')), ncol = 4, byrow = TRUE)
  cnv_gene_l_id <- matrix(unlist(strsplit(colnames(rt_cnv_gene_l), '[.]')), ncol = 4, byrow = TRUE)
  colnames(rt_cnv_gene_c) <- paste(cnv_gene_c_id[, 1], cnv_gene_c_id[, 2], cnv_gene_c_id[, 3], sep ='-')
  colnames(rt_cnv_gene_l) <- paste(cnv_gene_l_id[, 1], cnv_gene_l_id[, 2], cnv_gene_l_id[, 3], sep ='-')
  
  int_mut_cnv_cli <- intersect(intersect(row.names(rt_risk_cnv_order), colnames(rt_cnv_gene_c)), colnames(rt_cnv_gene_l))
  rt_risk_cnv_order_m <- rt_risk_cnv_order[match(int_mut_cnv_cli, row.names(rt_risk_cnv_order), nomatch = 0), ]
  rt_mvm_sur_all_omo <- rt_risk_cnv_order_m[order(rt_risk_cnv_order_m$risk_score), ]
  all(row.names(rt_mvm_sur_all_omo) == row.names(rt_risk_cnv_order_m))
  rt_cnv_gene_cm <- rt_cnv_gene_c[, match(int_mut_cnv_cli, colnames(rt_cnv_gene_c), nomatch = 0)]
  all(colnames(rt_cnv_gene_cm) == row.names(rt_mvm_sur_all_omo))
  rt_cnv_gene_lm <- rt_cnv_gene_l[, match(int_mut_cnv_cli, colnames(rt_cnv_gene_l), nomatch = 0)]
  all(colnames(rt_cnv_gene_lm) == row.names(rt_mvm_sur_all_omo))
  
  ###cnv plot
  rt_cnv_ttest <- as.matrix(rbind(rt_cnv_gene_lm, rt_mvm_sur_all_omo$risk_score))
  row.names(rt_cnv_ttest)[24777] <- 'risk_score'
  rt_cnv_ttest[c( which(rt_cnv_ttest == 0))] <- 'no_cnv'
  rt_cnv_ttest[c(which(rt_cnv_ttest == -1),  which(rt_cnv_ttest == 1), which(rt_cnv_ttest == -2), which(rt_cnv_ttest == 2))] <- 'cnv'
  rt_cnv_ttest <- data.frame(t(rt_cnv_ttest), stringsAsFactors = FALSE)
  rt_cnv_ttest$risk_score <- as.numeric(rt_cnv_ttest$risk_score)
  CountGeneMutRate(rt_cnv_ttest, mtype = 'cnv')
 

  CountPal <- function(gnam, rt_all){
    ###gnam with cnv
    ###colnames(rt_all) contain gene name and risk score
    
    rt_gene_score <-rt_all[, c(gnam, 'risk_score')]
    
    cnv_score <- as.numeric(rt_all$risk_score[grep("cnv", rt_all[, gnam])])
    no_cnv_score <- as.numeric(rt_all$risk_score[grep("no_cnv", rt_all[, gnam])])
    pval_list <-  t.test(cnv_score, no_cnv_score)
    pval = pval_list$p.value
    pval_nam <- c(gnam, pval)
    return(pval_nam)
  }
  
  gnam_list <- colnames(rt_cnv_ttest)[-24777]
  pvalue_list <- lapply(gnam_list, CountPal, rt_cnv_ttest)
  rt_cnv_p <- data.frame(do.call(rbind, pvalue_list), stringsAsFactors = FALSE)
  colnames(rt_cnv_p) <- c('gene_name', 'pval')
  rt_cnv_p_0.05 <- rt_cnv_p[which(rt_cnv_p$pval < 0.05), ]
  write.table(rt_cnv_p_0.05, file = 'DEG_genes_CNV.txt', row.names = FALSE, col.names = TRUE, sep = '\t')
  
  rt_cnv_gene_cm_p0.05 <- rt_cnv_gene_cm[match(rt_cnv_p_0.05$gene_name, row.names(rt_cnv_gene_cm), nomatch = 0), ]
  
  risk_score <- as.numeric(rt_cnv_ttest[, 'risk_score'])
  names(risk_score) <- row.names(rt_cnv_ttest)
  ha1 = HeatmapAnnotation(risk_score = anno_points(risk_score, gp = gpar(col = ifelse(risk_score > median(risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                          show_legend = TRUE,  annotation_height = unit( 30, "mm"), show_annotation_name = TRUE)
  
  pdf(file = paste(cancer, "cnv_gene_risk_score.pdf", sep = "_"), 7, 6)
  ht_list <- Heatmap(rt_cnv_gene_cm_p0.05, name = "copy number", col = colorRamp2(c(-2 , 0, 2), c(" blue", "white", "red")),
                     bottom_annotation = ha1, bottom_annotation_height = unit(3, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                     cluster_columns = FALSE, column_dend_reorder = FALSE, 
                     cluster_rows = TRUE, row_dend_reorder = FALSE, 
                     show_row_dend = TRUE, show_column_dend = FALSE,
                     show_row_names = FALSE, show_column_names = FALSE, row_names_side = "right") 
  draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()
  
}





