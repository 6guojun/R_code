RiskScoreHeat <- function(exp_mat, rt_cli_com, Fnam, col1 = col1, col2 = col2){
  ha1 = HeatmapAnnotation(risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), col1, col2)),width = unit(3, "cm"), axis = TRUE),
                          os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, col1, 'grey')),width = unit(3, "cm"), axis = TRUE),
                          
                          show_legend = c(TRUE, TRUE),
                          annotation_height = unit(c(30, 30), "mm"), show_annotation_name = TRUE)
  
  pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
  ht_list <- Heatmap(exp_mat, name = paste('expression level', 'z-score', sep = '\n'), col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(col2, "white", col1)),
                     bottom_annotation = ha1, bottom_annotation_height = unit(6, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                     cluster_columns = FALSE, column_dend_reorder = FALSE, 
                     cluster_rows = TRUE, row_dend_reorder = FALSE, 
                     show_row_dend = TRUE, show_column_dend = FALSE,
                     show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
  draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()
}
