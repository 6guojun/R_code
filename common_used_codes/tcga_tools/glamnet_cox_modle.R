##the function you can get the best cox model which is the most frequency 
##usage:  x_coef_value_list <- GLCoxMain(k, x, y, width, height, theme_coef)
##x is the expression table, the row.names is samples id and colname is the gene id
##y id the survival data which contain time and status
##k is the frequency times that you want to train in the model
##width is the boxplot width and height is the boxplot height
##theme_coef is theme of boxplot
##x_coef is expression matrix contained active covariates in the model 
##risk_genes with the active covariates
##coef_value_d is the coeffecient value of the risk genes
##shangjunv@163.com

library(ggplot2)
library(glmnet)
library(survival)

GLCoxMain <- function(k, x, y, width, height){
  ScreenKeyGS <- function(k){
    #x is the expression table, the row.names is samples id and colname is the gene id
    #y id the survival data which contain time and status
    #k is the frequency times that you want to train in the model
    #width is the boxplot width and height is the boxplot height
    #theme_coef is theme of boxplot
    #x_coef is expression matrix contained active covariates in the model 
    #risk_genes with the active covariates
    #coef_value_d is the coeffecient value of the risk genes
    
    cvfit = cv.glmnet(x, y, family = "cox")
    coef_min = coef(cvfit, s = "lambda.min")
    x_coef <- x[, which(coef_min[, 1] != 0)]
    risk_genes <- colnames(x_coef)
    coef_value_d <- coef_min[which(coef_min[, 1] != 0)]
    risk_genes_lam <- paste(risk_genes, collapse = ";")
    coef_list <- list(risk_genes_lam, x_coef, coef_value_d)
    return(coef_list)
  }
  
  CountGS <- function(t, risk_gens_vec, UG_list){
    pos_genes <- which(risk_gens_vec == UG_list[[t]])[1]
    leng_t <- length(which(risk_gens_vec == UG_list[[t]]))
    type_t <- paste(length(unlist(strsplit(UG_list[[t]], ';'))), 'genes', sep = '_')
    genes <- UG_list[[t]]
    GS_out <- c(type_t, genes, leng_t, pos_genes)
    return(GS_out)
  }
  
  coef_list <- sapply(1: k, ScreenKeyGS)
  risk_genes_list <- coef_list[seq(1, by = 3, length = k)]
  risk_gens_vec <- do.call(rbind, risk_genes_list)
  UG_list <- unique(risk_genes_list)
  
  t = 1: length(unique(risk_genes_list))
  GS_list <- lapply(t, CountGS, risk_gens_vec, UG_list)
  GS_mat <- data.frame(do.call(rbind, GS_list), stringsAsFactors = FALSE)
  colnames(GS_mat) <- c('type', 'genes', 'freuency', 'pos')
  GS_mat$type <- paste(GS_mat$pos, GS_mat$type, sep = '_')
  GS_mat$freuency <- as.numeric(GS_mat$freuency)
  GS_mat$pos <- as.numeric(GS_mat$pos)
  write.table(GS_mat, file = 'modle_gene_counts.txt', row.names = TRUE, col.names = TRUE, sep = "\t")
  
#get risk genes, risk genes expression table, risk genes' coefficient value
  
  GS_mat_d <- GS_mat[!(grep('genes', GS_mat$type) %in% grep('_0_genes', GS_mat$type)), ]
  print(GS_mat_d$freuency)
  risk_genes <- risk_genes_list[as.numeric(GS_mat_d$pos[unique(which(GS_mat_d$freuency == max(GS_mat_d$freuency)))])]
  x_coef_lsit <- coef_list[seq(2, by = 3, length = k)]
  x_coef <- x_coef_lsit[as.numeric(GS_mat_d$pos[which(GS_mat_d$freuency == max(GS_mat_d$freuency))])]
  coef_value_d_list <- coef_list[seq(3, by = 3, length = k)]
  coef_value_d <- coef_value_d_list[as.numeric(GS_mat_d$pos[which(GS_mat_d$freuency == max(GS_mat_d$freuency))])]
  
#do boxplot
  pdf(file = paste("fren_boxplot.pdf", sep = ""), width = width, height = height)
  p <- ggplot(GS_mat, aes(x = reorder(type, freuency) , y = freuency, fill = type)) + 
    geom_bar(stat = "identity")  + labs(title = 'the frequency of models') + ylab("frequency")  + xlab('type') +
    geom_text(aes(label = freuency),vjust = 1, hjust = 0.5, color = 'white', size = 8) 
  print(p)
  dev.off()
  
  x_coef_value_list <- list(risk_genes, x_coef, coef_value_d, GS_mat)
  return(x_coef_value_list)
}

