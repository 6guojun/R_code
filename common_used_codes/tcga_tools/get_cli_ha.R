###usage: GetCliHa(rt_cli_com, cancer)

GetCliHa <- function(rt_cli_com, cancer){
  if(cancer == "KIRP"){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(recurrence = rt_cli_com$recurrence, 
                            stage = rt_cli_com$stage, 
                            gender = rt_cli_com$gender,
                            col = list(recurrence= structure(names = c("1", "2", "3"), c('black', 'grey', 'white')),
                                       stage = structure(names = c("1", "2", "3"), c('grey', 'black', 'white')), 
                                       gender = structure(names = c("1", "2"), c( 'black', 'grey'))), 
                            
                            age = anno_points(rt_cli_com$age, gp = gpar(col = ifelse(rt_cli_com$age > median(rt_cli_com$age), "black", "grey")), width = unit(4, "cm"), axis = TRUE),
                            risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE),
                            os_time =  anno_points(rt_cli_com$OS_Time, gp = gpar(col = ifelse(rt_cli_com$OS_Status == 1, "black", "grey")),width = unit(4, "cm"), axis = TRUE),
                            
                            show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                            annotation_height = unit(c(5, 5, 5, 30, 30, 30), "mm"), show_annotation_name = TRUE)
    return(ha1)
  }
  if(cancer == "LAML"){
    print('element contain AML_FAB_subtype, AML_FAB_subtype, cytogenetic_risk_group, cytogenetic_abnormality, age, hemoglobin,
          blast_cell_BM, leukocyte, platelet, blast_cell_PB, risk_score, os_time')
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
  }
  if(cancer == "LUAD"){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE),
                            annotation_height = unit(30, "mm"), show_annotation_name = TRUE)
    return(ha1)
  }
  if(cancer == "LIHC"){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE),
                            annotation_height = unit(30, "mm"), show_annotation_name = TRUE)
    return(ha1)
  }
  if(cancer == "READ"){
    print('element contain recurrence, stage, gender, age, risk_score, os_time')
    print('age, risk_score, os_time must be numeric')
    ha1 = HeatmapAnnotation(risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE),
                            annotation_height = unit(30, "mm"), show_annotation_name = TRUE)
    return(ha1)
  } else {
      print('element contain recurrence, stage, gender, age, risk_score, os_time')
      print('age, risk_score, os_time must be numeric')
      ha1 = HeatmapAnnotation(risk_score = anno_points(rt_cli_com$risk_score, gp = gpar(col = ifelse(rt_cli_com$risk_score > median(rt_cli_com$risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE),
                              annotation_height = unit(30, "mm"), show_annotation_name = TRUE)
      return(ha1)
    }
}


 
    
    