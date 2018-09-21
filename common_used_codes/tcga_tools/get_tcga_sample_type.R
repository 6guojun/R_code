###the funciton was used to get tcga cancer or normal or all type and convert samples ID to patient ID
### #the format of sample id can be such as "TCGA.76.4927.01A", "TCGA-76-4927-01A" or "TCGA_76_4927_01A", you can choose based on split_type
### you must set sam_type to point which samples type you want to get 
###usage:
###rt_T <- split_tcga_tn(rt_sym, sam_type = "tumor")
###rt_N <- split_tcga_tn(rt_sym, sam_type = "normal")
###Email: shangjun@163.com
###20180412


###################################################
#function
#get tumor samples and make samples ID to patient ID
#get TCGA tumor or normal data
###################################################
split_tcga_tn <-function(rt_tcga, sam_type = c("tumor", "normal"), split_type = c("[.]", "-", "_")){
  #the colnames is samples id and row.name is gene id 
  #the format of sample id must be "TCGA.76.4927.01A" 
  #sam_type must be set
  get_sam_nam <- function(rt_tcga, split_type){
    if(split_type == "[.]"){
      sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "[.]")), ncol = 4, byrow = TRUE)
    }else if(split_type == "-"){
      sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "-")), ncol = 4, byrow = TRUE)
    } else if(split_type == "_") {
      sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "_")), ncol = 4, byrow = TRUE)
    } else {
      stop("you must set split_type ")
    }
    sam_nam <- paste(sam_G[, 1], sam_G[, 2], sam_G[, 3], sep = "-")
    sam_list <- list(sam_G, sam_nam)
    return(sam_list)
  }
 
   sam_list <- get_sam_nam(rt_tcga, "[.]")
   sam_G <- sam_list[[1]]
   sam_nam <- sam_list[[2]]
  
  if(sam_type == "tumor"){
    tcga_pos_T <- c(grep("01", sam_G[, 4]), grep("02", sam_G[, 4]), grep("03", sam_G[, 4]), grep("04", sam_G[, 4]), grep("05", sam_G[, 4]), grep("06", sam_G[, 4]), grep("07", sam_G[, 4]), grep("08", sam_G[, 4]), grep("09", sam_G[, 4]))
    rt_T <- rt_tcga[, tcga_pos_T]
    colnames(rt_T) <- sam_nam[tcga_pos_T]
    data <- rt_T
  } else if(sam_type == "normal"){
    tcga_pos_N <- c(grep("11", sam_G[, 4]))
    if(length(tcga_pos_N) > 0) {
      rt_N <- rt_tcga[, tcga_pos_N]
      colnames(rt_N) <- sam_nam[tcga_pos_N]
      data <- rt_N
    } else {
      warning("this is no normal samples in this cancer")
      next
    }
  } else {
    stop("you must set sam_type to point which type samples you want to get")
  }
}

####convert tcga samples to samples id and T/N group
GetSam <- function(mat_exp, ncol = 4, col_A = 'red', col_B = 'blue'){
  samples_id <- colnames(mat_exp)

  cancer_type <- matrix(unlist(strsplit(samples_id, '[.]')), ncol = ncol, byrow = TRUE)[, 4]
  cancer_type <- gsub("1..", "Normal", cancer_type)
  cancer_type <- gsub("0..", "Tumor", cancer_type)
  color_type <- gsub('Tumor', col_A, cancer_type)
  color_type <- gsub('Normal', col_B, color_type)
  rt_sam <- data.frame(cbind(samples_id, cancer_type, color_type), stringsAsFactors = FALSE)
  colnames(rt_sam) <- c('samples_id', 'group', 'color')
  return(rt_sam)
  
}
