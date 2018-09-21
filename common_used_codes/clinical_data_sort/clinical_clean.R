CliClean <- function(rt_cli, cancer){
  if(cancer == "LUAD"){
    rt_cli_T <- GetTumorS(rt_cli)
    rt_esc <- rt_cli_T[, c("submitter_id.samples", "age_at_initial_pathologic_diagnosis", "gender.demographic", "number_pack_years_smoked", "race.demographic",
                         "new_tumor_event_after_initial_treatment", "pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")]
    colnames(rt_esc) <- c("sample_id", "age", "gender", "smoke", "race", "recurrence", "T", "N", "M", "stage")
    
    ###convert clinical element to numeric
    rt_esc$gender[which(rt_esc$gender == "")] <- 'NA'
    rt_esc$race[which(rt_esc$race == "")] <- "not reported" 
    rt_esc$recurrence[which(rt_esc$recurrence == "")] <- "NA"
    rt_esc$T[grep("T1", rt_esc$T)] <- "T1"
    rt_esc$T[grep("T2", rt_esc$T)] <- "T2"
    rt_esc$T[grep("T3", rt_esc$T)] <- "T3"
    rt_esc$T[grep("T4", rt_esc$T)] <- "T4"
    rt_esc$T[which(rt_esc$T == "")] <- "TX"
    rt_esc$N[which(rt_esc$N == "")] <- "NX"
    rt_esc$M[which(rt_esc$M == "")] <- 'MX'
    rt_esc$stage[grep("stage iv", rt_esc$stage)] <- "stage IV"
    rt_esc$stage[grep("stage iii", rt_esc$stage)] <- "stage III"
    rt_esc$stage[grep("stage ii", rt_esc$stage)] <- "stage II"
    rt_esc$stage[grep("stage i", rt_esc$stage)] <- "stage I"
    rt_esc$stage[which(rt_esc$stage == "")] <- 'stage X'
  }
  rt_esc_c <- data.frame(apply(rt_esc, 2, as.character), stringsAsFactors = FALSE)
  colnames(rt_esc_c) <- colnames(rt_esc)
  write.table(rt_esc_c, file = paste(cancer, "_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  return(rt_esc_c)
}
