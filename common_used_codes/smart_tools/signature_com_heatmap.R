library(circlize)
library(ComplexHeatmap)

ComHeatSig <- function(mat_sur_score, risk_genes = risk_genes, risk_title = risk_title){
  
  risk_order <- order(mat_sur_score[, risk_title])
  risk_score <-mat_sur_score[, risk_title][risk_order]
  os_time <- mat_sur_score$OS_Time[risk_order]
  status <- mat_sur_score$OS_Status[risk_order]
  
  exp <- mat_sur_score[, risk_genes][risk_order, ]
  exp_scale <- t(scale(exp))
  
  pdf(file = paste(risk_title, "_com_cluster.pdf", sep = ""), 10, 8)
  ha1 = HeatmapAnnotation(risk_score = anno_points(as.numeric(risk_score), gp = gpar(col = ifelse( risk_score > median(risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE), 
                          os_time =  anno_points(as.numeric(os_time), gp = gpar(col = ifelse( status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                          show_legend = c(TRUE, TRUE),
                          annotation_height = unit(c(30, 30), "mm"), show_annotation_name = TRUE)
  
  ht_list <- Heatmap(exp_scale, name = "center scaled expression", col = colorRamp2(c(min(exp_scale)/2 , median(exp_scale), max(exp_scale)/2), c(" blue", "white", "red")),
                     bottom_annotation = ha1, bottom_annotation_height = unit(6, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                     cluster_columns = FALSE, column_dend_reorder = FALSE, 
                     cluster_rows = TRUE, row_dend_reorder = FALSE, 
                     show_row_dend = TRUE, show_column_dend = FALSE,
                     show_row_names = TRUE, show_column_names = FALSE) 
  
  draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()
}
