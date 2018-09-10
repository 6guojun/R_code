library(ComplexHeatmap)
library(circlize)

expr = readRDS(paste0(system.file(package = "ComplexHeatmap"), "/extdata/gene_expression.rds"))
mat = as.matrix(expr[, grep("cell", colnames(expr))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))

type = gsub("s\\d+_", "", colnames(mat))
age = c(1:24)
ha1 = HeatmapAnnotation(age = anno_points(age, gp = gpar(col = ifelse(age > 5, "black", "red")), axis = TRUE), 
                                         annotation_height = unit(c(30, 5), "mm"))
ha2 = HeatmapAnnotation(type = type)

ht_list <- Heatmap(mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                   top_annotation = ha2, bottom_annotation = ha1, top_annotation_height = unit(4, "mm"), 
        clustering_distance_rows = "euclidean", row_dend_reorder = TRUE,
        clustering_distance_columns = "spearman", column_dend_reorder = TRUE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        show_row_names = FALSE, show_column_names = FALSE) +
  Heatmap(base_mean, name = "base_mean", show_row_names = FALSE, width = unit(5, "mm")) +
  Heatmap(expr$length, name = "length", col = colorRamp2(c(0, 1000000), c("white", "orange")),
          heatmap_legend_param = list(at = c(0, 200000, 400000, 60000, 800000, 1000000), 
                                      labels = c("0kb", "200kb", "400kb", "600kb", "800kb", "1mb")),
          width = unit(5, "mm")) +
  Heatmap(expr$type, name = "type", width = unit(5, "mm"))

draw(ht_list, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
