###convet a matrix with expression  to matrix with gene pairs 

script_dir='/Users/stead/Documents/SourceTree/R/common_used_codes/'
source(paste(script_dir, 'tcga_tools/get_exp_sur_cli_merge_data.R', sep = ""))
source(paste(script_dir, 'tcga_tools/get_tcga_sample_type.R', sep = ""))
source(paste(script_dir, 'smart_tools/EA2SM.R', sep = ""))


###convert expression matrix to gene pairs matrix which contain only 1 or 0
Exp2GPS <- function(mat_exp, rm_cons = c('yes', 'no')){
  ###prepare a data matrxi which contain rownames (gene id) and colnames (samples id)
  GetPairMat <- function(num_i, mat_exp, gene_comb){
    mat_d <-  data.frame(mat_exp[gene_comb[num_i, 1], ] - mat_exp[gene_comb[num_i, 2], ], stringsAsFactors = FALSE)
    row.names(mat_d) <- paste(gene_comb[num_i, 1], gene_comb[num_i, 2], sep = '_')
    return(mat_d)
  }
  ###get gene combn
  gene_list <- row.names(mat_exp)
  gene_comb_p <- t(combn(gene_list, 2))
  num_i <- 1: length(gene_comb_p[, 1])
  mat_pair_list <- lapply(num_i, GetPairMat, mat_exp, gene_comb_p)
  mat_pair <- do.call(rbind, mat_pair_list) 
  if(rm_cons == 'yes'){
    ###sort data
    mat_pair_d <- mat_pair[!apply(mat_pair, 1, function(x){all(x > 0)| all(x < 0) | all(x = 0)}), ]
    mat_pair_d[mat_pair_d >= 0] <- 0
    mat_pair_d[mat_pair_d < 0] <- 1
    return(mat_pair_d)
  } else if (rm_cons == 'no'){
    mat_pair[mat_pair >= 0] <- 0
    mat_pair[mat_pair < 0] <- 1
    return(mat_pair)
  }
}

###get gene pair matrix with survival data
PreGPS <- function(data_dir = data_dir, cancer_type = cancer_type, dataset_id = c('TCGA', "GEO"), rm_cons = c('yes', 'no'), tar_genes = tar_genes, GSE_ID = GSE_ID, 
                   in_type = in_type, out_type = out_type, dup_type = dup_type) {
  ###TCGA dataset
  if(dataset_id == 'TCGA'){
    rt_exp <- read.table(file = paste(data_dir, '/', cancer_type, '/RNA_seq/TCGA-', cancer_type, '.htseq_fpkm.tsv', sep = ''), header = TRUE, 
                         row.names = 1, sep = '\t', stringsAsFactors = FALSE)
    rt_sur <- read.table(file = paste(data_dir, '/', cancer_type, "/phenotype/TCGA-", cancer_type, ".survival.tsv", sep = ""), header = TRUE,
                         sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
    
    rt_exp_sym <- EA2SM(rt_exp, in_type = in_type, out_type = out_type, dup_type = dup_type)  #conver ensemble id to gene name
    rt_exp_int <- rt_exp_sym[match(tar_genes, row.names(rt_exp_sym), nomatch = 0), ] #get target genes expression table
    
    #convert expressino table to gene pairs table
    if(rm_cons == 'yes'){
      mat_pair <- Exp2GPS(rt_exp_int, rm_cons = 'yes')
    } else if (rm_cons == 'no'){
      mat_pair <- Exp2GPS(rt_exp_int, rm_cons = 'no')
    }
    ###get exp and sur merge data
    mat_pair_T <- split_tcga_tn(mat_pair, sam_type = "tumor")   #achive tumor sample
    exp_sur_list <- TCGAExpSurM(mat_pair_T, rt_sur)
    mat_pair_m <- exp_sur_list[[1]]
    mat_sur_m <- exp_sur_list[[2]]
    pair_nam <- row.names(mat_pair_m)
    mat_pair_sur <- data.frame(cbind(t(mat_pair_m), mat_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)
    colnames(mat_pair_sur) <- c(pair_nam, "OS_Status", "OS_Time")
    return(mat_pair_sur)
    
    ###GEO dataset###
  } else if (dataset_id == 'GEO') {
    rt_exp <- read.table(file = paste(data_dir, '/', cancer_type, '/', GSE_ID, '/expression_data/', GSE_ID, '_mat.txt', sep = '')
                         , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rt_sur <- read.table(file = paste(data_dir, '/', cancer_type, '/', GSE_ID, '/survival_data/survival_data.txt', sep = ''), 
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    rt_exp_sym <- EA2SM(rt_exp, in_type = in_type, out_type = out_type, dup_type = dup_type)   #conver ensemble id to gene name
    rt_exp_int <- rt_exp_sym[match(tar_genes, row.names(rt_exp_sym), nomatch = 0), ] #get target genes expression table
    #convert expressino table to gene pairs table
    if(rm_cons == 'yes'){
      mat_pair <- Exp2GPS(rt_exp_int, rm_cons = 'yes')
    } else if (rm_cons == 'no'){
      mat_pair <- Exp2GPS(rt_exp_int, rm_cons = 'no')
    }
    exp_sur_list <- GEOExpCliM(mat_pair, rt_sur)
    mat_pair_m <- exp_sur_list[[1]]
    mat_sur_m <- exp_sur_list[[2]]
    pair_nam <- row.names(mat_pair_m)
    mat_pair_sur <- data.frame(cbind(t(mat_pair_m), mat_sur_m[, c("OS_Status", "OS_Time")]), stringsAsFactors = FALSE)
    colnames(mat_pair_sur) <- c(pair_nam, "OS_Status", "OS_Time")
    return(mat_pair_sur)
    
  } else {
    stop('dataset_id should contain TCGA or GSE' )
  }
}

GetGS <- function(data_dir = data_dir, cancer_type = cancer_type, dataset_id = c('TCGA', "GEO"), tar_genes = tar_genes, GSE_ID = GSE_ID, 
                   in_type = in_type, out_type = out_type, dup_type = dup_type) {
  ###TCGA dataset
  if(dataset_id == 'TCGA'){
    rt_exp <- read.table(file = paste(data_dir, '/', cancer_type, '/RNA_seq/TCGA-', cancer_type, '.htseq_fpkm.tsv', sep = ''), header = TRUE, 
                         row.names = 1, sep = '\t', stringsAsFactors = FALSE)
    rt_sur <- read.table(file = paste(data_dir, '/', cancer_type, "/phenotype/TCGA-", cancer_type, ".survival.tsv", sep = ""), header = TRUE,
                         sep = "\t", row.names = NULL, stringsAsFactors = FALSE)
    
    rt_exp_sym <- EA2SM(rt_exp, in_type = in_type, out_type = out_type, dup_type = dup_type)  #conver ensemble id to gene name
    rt_exp_int <- rt_exp_sym[match(tar_genes, row.names(rt_exp_sym), nomatch = 0), ] #get target genes expression table

    ###get exp and sur merge data
    mat_pair_T <- split_tcga_tn(rt_exp_int, sam_type = "tumor")   #achive tumor sample
    exp_sur_list <- TCGAExpSurM(mat_pair_T, rt_sur)
    mat_pair_m <- exp_sur_list[[1]]
    mat_sur_m <- exp_sur_list[[2]]
    pair_nam <- row.names(mat_pair_m)
    mat_pair_sur <- data.frame(cbind(t(mat_pair_m), mat_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)
    colnames(mat_pair_sur) <- c(pair_nam, "OS_Status", "OS_Time")
    return(mat_pair_sur)
    
    ###GEO dataset###
  } else if (dataset_id == 'GEO') {
    rt_exp <- read.table(file = paste(data_dir, '/', cancer_type, '/', GSE_ID, '/expression_data/', GSE_ID, '_mat_log.txt', sep = '')
                         , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    rt_sur <- read.table(file = paste(data_dir, '/', cancer_type, '/', GSE_ID, '/survival_data/survival_data.txt', sep = ''), 
                         header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    rt_exp_sym <- EA2SM(rt_exp, in_type = in_type, out_type = out_type, dup_type = dup_type)   #conver ensemble id to gene name
    rt_exp_int <- rt_exp_sym[match(tar_genes, row.names(rt_exp_sym), nomatch = 0), ] #get target genes expression table
  
    exp_sur_list <- GEOExpCliM(rt_exp_int, rt_sur)
    mat_pair_m <- exp_sur_list[[1]]
    mat_sur_m <- exp_sur_list[[2]]
    pair_nam <- row.names(mat_pair_m)
    mat_pair_sur <- data.frame(cbind(t(mat_pair_m), mat_sur_m[, c("OS_Status", "OS_Time")]), stringsAsFactors = FALSE)
    colnames(mat_pair_sur) <- c(pair_nam, "OS_Status", "OS_Time")
    return(mat_pair_sur)
    
  } else {
    stop('dataset_id should contain TCGA or GSE' )
  }
}


