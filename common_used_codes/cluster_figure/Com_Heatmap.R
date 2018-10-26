#do heatmap
SubHeat <- function(group_sub, rt_sub_ord, data_type, method_sub, num_k, color_type, 
                    cluster_columns = FALSE, show_column_dend = FALSE,  column_dend_reorder = FALSE, show_column_names = FALSE, 
                    cluster_rows = TRUE, show_row_dend = TRUE, row_dend_reorder = TRUE, show_row_names = FALSE){
  
  print('group_sub and the expression mat must have the same samples names')
  ha1 = HeatmapAnnotation(group_sub = group_sub, 
                          col = list(group_sub= structure(names = unique(group_sub), color_type[1:length(unique(group_sub))])))
  pdf(file = paste(data_type, method_sub, num_k, 'ht.pdf', sep = '_'), 10, 8)
  ht_list <- Heatmap(rt_sub_ord, name = "expression level \n z-score",
                     col = colorRamp2(c(min(rt_sub_ord)/3 , (min(rt_sub_ord) + max(rt_sub_ord))/3, max(rt_sub_ord)/2), c(" blue", "white", "red")),
                     top_annotation = ha1, top_annotation_height = unit(0.5, "cm"), 
                     cluster_columns = cluster_columns, column_dend_reorder = column_dend_reorder, 
                     cluster_rows = cluster_rows, row_dend_reorder = row_dend_reorder, 
                     show_row_dend = show_row_dend, show_column_dend = show_column_dend,
                     show_row_names = show_row_names, show_column_names = show_column_names, row_names_side = "right") 
  draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()
  
  print('heatmap finished')
}
