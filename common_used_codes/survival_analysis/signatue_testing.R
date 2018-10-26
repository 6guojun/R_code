#the colnames of mat_testing_sur contain genes and survival information
#the rownames of mat_testing_sur contain samples id
#risk genes is the signature genes
#coef_value is the coef value which count through cox
#dataset_id is used to set file names
#count GSE37745 risk score

script_dir <- '/Users/stead/Documents/SourceTree/R/common_used_codes/' 
source(paste(script_dir, 'tcga_tools/time_ROC.R', sep = ""))
source(paste(script_dir, 'smart_tools/signature_com_heatmap.R', sep = ""))
source(paste(script_dir, 'survival_analysis/survival_analysis.R', sep = ""))

SignatureTest <- function(mat_testing_sur = mat_testing_sur, risk_genes = risk_genes, coef_value = coef_value, dataset_id = dataset_id, 
                          feature = c('OS', 'DFS'), exp_type = exp_type){
  #
  mat_exp_testing <- mat_testing_sur[, risk_genes]
  risk_score_testing <- apply(t(mat_exp_testing)*coef_value, 2, function(x){sum(x)})
  mat_testing_sur$risk_score <- risk_score_testing
  if(feature == 'OS'){
  mat_testing_score <- mat_testing_sur[, c(risk_genes, 'risk_score', 'OS_Status', 'OS_Time')]
  title_name <- paste(dataset_id, '_risk_score', sep = "")
  colnames(mat_testing_score) <- c(risk_genes, title_name, 'OS_Status', 'OS_Time')
  } else if (feature == 'DFS'){
  mat_testing_score <- mat_testing_sur[, c(risk_genes, 'risk_score', 'DFS_Status', 'DFS_Time')]
  title_name <- paste(dataset_id, '_risk_score', sep = "")
  colnames(mat_testing_score) <- c(risk_genes, title_name, 'DFS_Status', 'DFS_Time')
  }
  lapply(title_name, SurMainFun, mat_testing_score, DatType = "ConType", feature = feature, OutType = "SurF", SurvType = 'ALL', len_a = 'high_risk', len_b = 'low_risk')
  
  ###time dependent ROC
  if(feature == 'OS'){
  CountTimeROC(mat_testing_score, 'OS_Time', 'OS_Status', title_name)
  } else if (feature == 'DFS'){ 
  CountTimeROC(mat_testing_score, 'DFS_Time', 'DFS_Status', title_name)  
  }
  
  ###complex heatmap
  ComHeatSig(mat_testing_score, risk_genes = risk_genes, risk_title = title_name, exp_type = exp_type)
  return(mat_testing_score)
}
