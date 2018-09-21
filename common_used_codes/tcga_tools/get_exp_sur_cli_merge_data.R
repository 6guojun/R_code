###you will get a list which contain expression table, clinical data and survival data 
###all of them contain the same sample id and the same order
###usege: list_luad <- get_exp_sur_cli(rt, rt_cli, rt_sur)


GetExpSurCli <- function(rt_T, rt_cli, rt_sur){
  print('the colnames of expression table should contain "Ensembl_ID", 
         you will get a list containing expression table, clinical data and survival data
         all of which contain the same sample id and the same order, 
         the expression data contain only tumor samples')
  
  ###survival data
  sam_sur <- matrix(unlist(strsplit(rt_sur$sample, "-")), ncol = 4, byrow = TRUE)
  rt_sur$sample <- paste(sam_sur[, 1], sam_sur[, 2], sam_sur[, 3], sep = "-")
  
  #clinical data
  print("the first column must be samples id which was combined wiht '-'")
  sam_cli <- matrix(unlist(strsplit(rt_cli[, 1], "-")), ncol = 4, byrow = TRUE)
  rt_cli$submitter_id.samples <- paste(sam_cli[, 1], sam_cli[, 2], sam_cli[, 3], sep = "-")
  
  ###intersection sampples
  sam_m <- intersect(colnames(rt_T),intersect(rt_sur$sample, rt_cli$submitter_id.samples))
  
  pair_T <- match(sam_m, colnames(rt_T), nomatch = NA)
  pair_S <- match(sam_m, rt_sur$sample , nomatch = NA)
  pair_C <- match(sam_m, rt_cli$submitter_id.samples , nomatch = NA)
  
  rt_T_m <- rt_T[, pair_T]
  rt_S_m <- rt_sur[pair_S, ]
  rt_C_m <- rt_cli[pair_C, ]
  
  esc_list <- list(rt_T_m,  rt_C_m, rt_S_m)
  if(all(colnames(rt_T_m) == rt_S_m$sample)&all(rt_S_m$sample == rt_C_m$submitter_id.samples)){
    return(esc_list)
  } else {
    stop("the sample id is diff among expression data, survival data and clinical data ")
  }
}


GEOExpCliM <- function(rt_exp, rt_cli_sur){
  int_id <- intersect(colnames(rt_exp), rt_cli_sur$Accession)
  rt_exp_m <- rt_exp[, match(int_id, colnames(rt_exp), nomatch = 0)]
  rt_cli_sur_m <- rt_cli_sur[match(int_id, rt_cli_sur$Accession, nomatch = 0), ]
  esc_list <- list(rt_exp_m,  rt_cli_sur_m)
  return(esc_list)
}


TCGAExpSurM <- function(rt_exp, rt_sur){
  
  
  ###survival data
  sam_sur <- matrix(unlist(strsplit(rt_sur$sample, "-")), ncol = 4, byrow = TRUE)
  rt_sur$sample <- paste(sam_sur[, 1], sam_sur[, 2], sam_sur[, 3], sep = "-")

  ###intersection sampples
  int_id <- intersect(colnames(rt_exp), rt_sur$sample)
  
  rt_exp_m <- rt_exp[, match(int_id, colnames(rt_exp), nomatch = 0)]
  rt_sur_m <- rt_sur[match(int_id, rt_sur$sample , nomatch = 0), ]

  esc_list <- list(rt_exp_m,  rt_sur_m)
  if(all(colnames(rt_exp_m) == rt_sur_m$sample)){
    return(esc_list)
  } else {
    stop("the sample id is diff among expression data, survival data and clinical data ")
  }
}

MergeCliOth <- function(rt_cli, rt_other){
  sam_int <- intersect(rt_cli$sample_id, row.names(rt_other))
  rt_cli_m <- rt_cli[match(sam_int, rt_cli$sample_id, nomatch = 0), ]
  rt_other_m <- rt_other[match(sam_int, row.names(rt_other), nomatch = 0), ]
  rt_other_cli <- data.frame(cbind(rt_other_m, rt_cli_m), stringsAsFactors = FALSE)
  return(rt_other_cli)
}

