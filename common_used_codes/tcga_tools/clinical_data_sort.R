CliSort <- function(rt_cli, cancer){
  if(cancer == "LUAD"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", "number_pack_years_smoked", "race.demographic",
                         "new_neoplasm_event_type", "pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "smoke", "race", "recurrence", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))#female male
    rt_esc$race <- factor(rt_esc$race, labels=c(1:5))#american indian or alaska native asian black or african american not reported white
    rt_esc$recurrence <- factor(rt_cli_m$new_tumor_event_after_initial_treatment, labels=c(3, 2, 1))#''  NO YES
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 4, 5))#T1 T1a T1b T2 T2a T2b T3 T4 TX
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 3))#''  N0 N1 N2 N3 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2, 3))#''  M0 M1 M1a M1b MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2))#not reported stage i stage ia stage ib stage ii stage iia stage iib stage iiia stage iiib stage iv
  }
  if(cancer == "UVM"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses",  "new_tumor_event_after_initial_treatment")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", "recurrence")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3))
    rt_esc$N <- factor(rt_esc$N, labels=c(2, 1, 2))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2))
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 1, 2))
  }
  if(cancer == "BRCA"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5))
    rt_esc$N <- factor(rt_esc$N, labels=c(rep(1, 4), rep(2, 11), 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(1, 2, 3, 1))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3))
  }
  if(cancer == "ESCA"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id",  "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(6, 1, 2, 3, 4, 5, 5))
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
  }
  if(cancer == "HNSC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(6, 1, 2, 3, 4, 5, 5, 5, 6))
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 2, 2, 2, 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2))
  }
  if(cancer == "KIRP"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", "new_tumor_event_after_initial_treatment")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", "recurrence")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))# female male
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5))# T1 T1a T1b T2 T2a T2b T3 T3a T3b T3c T4 TX
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 3))#"" N0 N1 N2 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))#"" M0 M1 MX
    rt_esc$stage <- factor(rt_esc$stage, labels = c(3, 1, 1, 2, 2))# not reported stage i stage ii stage iii stage iv
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels = c(3, 2, 1))#  NO YES
  }
  if(cancer == "KIRC"){
    rt_esc <-  rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4))
    rt_esc$N <- factor(rt_esc$N, labels=c(1, 2, 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
    rt_esc$stage <- factor(rt_esc$stage, labels = c(3, 1, 1, 2, 2))
  }
  if(cancer == "LAML"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_diagnosis.diagnoses", "gender.demographic", "lab_procedure_hemoglobin_result_specified_value",
                         "lab_procedure_bone_marrow_blast_cell_outcome_percent_value", "lab_procedure_leukocyte_result_unspecified_value", "platelet_result_count", 
                         "lab_procedure_blast_cell_outcome_percentage_value", "leukemia_french_american_british_morphology_code", 
                         "acute_myeloid_leukemia_calgb_cytogenetics_risk_category", "cytogenetic_abnormality")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "hemoglobin", "blast_cell_BM", "leukocyte", "platelet", "blast_cell_PB", 
                          "AML_FAB_subtype", "cytogenetic_risk_group", "cytogenetic_abnormality")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels = c(1, 2))#female male
    rt_esc$AML_FAB_subtype <- factor(rt_esc$AML_FAB_subtype, labels = c(1, 2, 3, 4, 5, 6, 7, 8)) #"M0 Undifferentiated" M1 M2 M3 M4 M5 M6 M7 Not Classified
    rt_esc$cytogenetic_risk_group <- factor(rt_esc$cytogenetic_risk_group, labels = c(3, 1, 2, 2)) #"" "Favorable" "Intermediate/Normal" "Poor"
    rt_esc$cytogenetic_abnormality <- factor(rt_esc$cytogenetic_abnormality, labels = c(3, 2, 2, 1, 2, 2, 2, 2, 2, 2)) #" " "+8" "Complex - >= 3 distinct abnormalities" "Normal" "del(5q) / 5q-" "del(7q) / 7q-" "inv(16)" "t(15;17) and  variants" "t(8;21)" "t(9;11)"
    return(rt_esc)
  } 
  if(cancer == "LGG"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", "neoplasm_histologic_grade")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "grade")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$grade <- factor(rt_esc$grade, labels=c(3, 1, 2))
  } 
  if(cancer == "OV"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "clinical_stage")]
    colnames(rt_esc) <- c("sample_id", "age", "clinical_stage")
    
    ###convert clinical element to numeric
    rt_esc$clinical_stage <- factor(rt_esc$clinical_stage, labels=c(3, 1, 1, 1, 1, 2, 2, 2, 2))
  } 
  if(cancer == "SKCM"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
    rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 2, 2, 3, 3, 4, 4, 4))
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, rep(2, 8), 3))
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2))
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2))
  } 
  if(cancer == "LUSC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 4)) # T1 T1a T1b T2 T2a T2b T3 T4
    rt_esc$N <- factor(rt_esc$N, labels=c(1, rep(2, 3), 3)) #N0 N1 N2 N3 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2, 3)) #  M0 M1 M1a M1b MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2)) # not reported stage i stage ia stage ib stage ii stage iia stage iib stage iii stage iiia stage iiib stage iv
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #  NO YES
  } 
  if(cancer == "LIHC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 2, 2, 3, 3, 3, 4, 5)) #'' T1 T2 T2a T2b T3 T3a T3b T4 TX
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 3)) #'' N0 N1 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(1, 2, 3)) # M0 M1 MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2, 2, 2, 2)) #  not reported stage i stage ii stage iii stage iiia stage iiib stage iiic stage iv stage iva stage ivb
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  if(cancer == "READ"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 3, 4, 4, 4)) #'' T1 T2 T3 T4 T4a T4b
    rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, rep(2, 7), 3)) #'' N0 N1 N1a N1b N1c N2 N2a N2b NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3)) #'' M0 M1 M1a MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, rep(1, 5), rep(2, 6))) #not reported stage i stage ii stage iia stage iib stage iic stage iii stage iiia stage iiib stage iiic stage iv stage iva
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  if(cancer == "SARC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", 'new_tumor_event_after_initial_treatment', 'tumor_depth')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", 'recurrence', 'tumor_depth')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$tumor_depth <- factor(rt_cli_m$tumor_depth, labels=c(3, 2, 1))#'' Deep Superficial
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  if(cancer == "KIRC"){
    rt_esc <- rt_cli[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", 
                         "pathologic_N", "pathologic_M", "tumor_stage.diagnoses", 'new_tumor_event_after_initial_treatment')]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "T", "N", "M", "stage", 'recurrence')
    
    ###convert clinical element to numeric
    rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2)) #female male
    rt_esc$T <- factor(rt_esc$T, labels=c(5, rep(1, 2), rep(2, 3), rep(3, 4), 4)) #'' T1 T1a T1b T2 T2a T2b T3 T3a T3b T3c T4
    rt_esc$N <- factor(rt_esc$N, labels=c(1, 2, 3)) # N0 N1 NX
    rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3)) #'' M0 M1 MX
    rt_esc$stage <- factor(rt_esc$stage, labels=c(3, rep(1, 2), rep(2, 2))) # not reported stage i stage ii stage iii stage iv
    rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 2, 1)) #'' NO YES
  }
  rt_esc_c <- data.frame(apply(rt_esc, 2, as.character), stringsAsFactors = FALSE)
  colnames(rt_esc_c) <- colnames(rt_esc)
  write.table(rt_esc_c, file = paste(cancer, "_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  return(rt_esc_c)
}


