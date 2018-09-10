###get univaiate analysis and LogRank gene
LogRankUcoxSigGs <- function(mat_exp_sur, data_id, DatType= c('ConType', 'LogType'), script_dir){
  print('change the content of survival_analysis_2.R')
  print('mat_exp_sur contain OS_Time and OS_Status in the end of mat_exp_sur colnames')
  source(paste(script_dir, 'survival_analysis/survival_analysis.R', sep = ""))
  unam <- colnames(mat_exp_sur)[1: (dim(mat_exp_sur)[2]-2)]
  uvm_list <- lapply(unam, uvm_count, mat_exp_sur)
  data_uvm <- do.call(rbind, uvm_list)
  data_uvm <- na.omit(data_uvm)
  sig_GS <-  row.names(data_uvm)[data_uvm$pvalue < 0.05]
  rt_uvm_sig <- mat_exp_sur[, match(sig_GS, colnames(mat_exp_sur), nomatch = 0)]
  
  ###do survival analysis
  gnam_log <- colnames(rt_uvm_sig)
  sur_plist <- lapply(gnam_log, SurMainFun, mat_exp_sur, DatType = DatType, feature = "OS", OutType = "SurP", SurvType = 'ALL')
  data_sur <- do.call(rbind, sur_plist)
  sig_GS_m <- data_sur[, 1][which(as.numeric(data_sur[, 2]) < 0.05 )]
  data_uvm_logRank_sig <- data_uvm[sig_GS_m, ]
  write.table(data_uvm_logRank_sig, file = paste(data_id, '_uvm_logRank_sig.txt', sep = ''), row.names = TRUE, col.names = TRUE, sep = '\t')
  return(sig_GS_m)
}
