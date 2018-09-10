###your rt_cli_com data should contain 
###do complex heatmap
PlotCliComHeat <- function(exp_mat, rt_cli_com, cancer, Fnam){
  if(cancer == "KIRP"){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'LAML'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(gender = rt_cli_com$gender, AML_FAB_subtype = rt_cli_com$AML_FAB_subtype, cytogenetic_risk_group = rt_cli_com$cytogenetic_risk_group,
                            cytogenetic_abnormality = rt_cli_com$cytogenetic_abnormality, hemoglobin = rt_cli_com$hemoglobin,  blast_cell_BM = rt_cli_com$blast_cell_BM,
                            leukocyte = rt_cli_com$leukocyte, platelet = rt_cli_com$platelet, blast_cell_PB = rt_cli_com$blast_cell_PB, 
                            col = list(gender = structure(names = c("1", "2"), c( 'black', 'grey')), 
                                       AML_FAB_subtype = structure(names = c("1", "2", "3", "4", "5", "6", "7", "8"), c("green4", "gold3", "orange1", "firebrick2", "plum3", "gold2", "darkred",  "yellow2")), 
                                       cytogenetic_risk_group = structure(names = c("1", "2", "3"), c("grey", "black", "white")),
                                       cytogenetic_abnormality = structure(names = c("1", "2", "3"), c("grey", "black", "white")), 
                                       hemoglobin = colorRamp2(c(0, max(rt_cli_com$hemoglobin)), c("white", "red")), 
                                       blast_cell_BM = colorRamp2(c(0, max(rt_cli_com$blast_cell_BM)), c("white", "red")), 
                                       leukocyte = colorRamp2(c(0, max(rt_cli_com$leukocyte)), c("white", "red")), 
                                       platelet = colorRamp2(c(0, max(rt_cli_com$platelet)), c("white", "red")), 
                                       blast_cell_PB = colorRamp2(c(0, max(rt_cli_com$blast_cell_PB)), c("white", "red"))),
    
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 5, 5, 5, 5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 10)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(13.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'LUAD'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'LIHC'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(2.5, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(2.5, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(2.5, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 25, 25, 25), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 9)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)/1.2), c(" blue", "white", "firebrick2")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(9, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
    
  }
  if(cancer == 'READ'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'SARC'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            tumor_depth = rt_cli_com$tumor_depth, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       tumor_depth = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
  if(cancer == 'KIRC'){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > 60, "black", "grey")), width = unit(3, "cm"),  lty = 3, axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(3, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(3, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    
    pdf(file = paste(cancer, Fnam, "risk_score.pdf", sep = "_"), 10, 8)
    ht_list <- Heatmap(exp_mat, name = "center scaled expression", col = colorRamp2(c(min(exp_mat) , 0, max(exp_mat)), c(" blue", "white", "red")),
                       bottom_annotation = ha1, bottom_annotation_height = unit(10.5, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                       cluster_columns = FALSE, column_dend_reorder = FALSE, 
                       cluster_rows = TRUE, row_dend_reorder = FALSE, 
                       show_row_dend = TRUE, show_column_dend = FALSE,
                       show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right") 
    draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
    dev.off()
  }
}
  