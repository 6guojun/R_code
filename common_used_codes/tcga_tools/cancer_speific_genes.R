### screenig cancer specific genes with t-test, pvalue < 0.05 and log2FC > 1
### GP1 and GP2 are two cancer name
### SGs_list <- SGTest(GP1, GP2, rt, rt_sam, FC, var.equal = var.equal, paired = FALSE)

SGTest <- function(GP1, GP2, rt, rt_sam, FC, var.equal = var.equal, paired = FALSE){
  #you can get one table which contain information of pvalue and fold change
  #you can aslo get one expression table of DEG 
  #you can choice achive two table as above
  print('GP1 must the cancer contained specific genes')
  print('the length and order of rt colnames and samples from rt_sam must be the same')
  
  DrawMat <- function(GP1, GP2, rt, rt_sam){
    #GP2 and GP2 are the name of your group
    #rt is the expression table 
    #rt_sam contain samples name and its group name
    #this function will get a new expression table which show the group information in colnames
    
    colnames(rt) <- paste(colnames(rt), rt_sam$group, sep = "_")
    GP1_N <- grep(GP1, colnames(rt))
    GP2_N <- grep(GP2, colnames(rt))
    data_T <- data.frame(rt[, c(GP1_N, GP2_N)], stringsAsFactors = FALSE)
    return(data_T)
  }
  
  
  GetPval <- function(Gnam, GP1, GP2, data, var.equal = var.equal, paired = FALSE){
    #get the pvalue one gene of T-test  between two group
    
    data_Gnam <- data[Gnam, ]
    GPx <- c(grep(GP1, colnames(data_Gnam)))
    GPy <- c(grep(GP2, colnames(data_Gnam)))
    x <- as.numeric(data_Gnam[, GPx])
    y <- as.numeric(data_Gnam[, GPy])
    log2FC <- mean(x) - mean(y)
    if(all(x[1] == (c(x, y)))){
      return(NA)
    } else{
      Obj <- t.test(x, y, var.equal = var.equal, paired = FALSE)
      Pval <- Obj$p.value
      Padj <-  p.adjust(Pval, method = "fdr")
      Cval <- c(Pval, Padj, log2FC)
      return(Cval) 
    }
  }
  
  
  MainTest <- function(GP1, GP2, rt, rt_sam, var.equal = var.equal, paired = FALSE){
    #get one table which contain information of pvalue and fold change
    
    rt_data <- DrawMat(GP1, GP2, rt, rt_sam)
    Gnam <- row.names(rt_data)
    Cval_L <- lapply(Gnam, GetPval, GP1, GP2, rt_data, var.equal = var.equal, paired = FALSE)
    data_Pval_FC <- data.frame(do.call(rbind, Cval_L), stringsAsFactors = FALSE)
    row.names(data_Pval_FC) <- row.names(rt_data)
    colnames(data_Pval_FC) <- c("Pval", "Padj", "log2FC")
    data_Pval_FC$gene_id <- row.names(data_Pval_FC)
    data_Pval_FC$group <- paste(GP1, GP2, sep = "_")
    return(data_Pval_FC)
  }
  
  ###get the Pval and FC table
  exp_mat <- DrawMat(GP1, GP2, rt, rt_sam)
  aLL_Pval_FC <- MainTest(GP1, GP2, rt, rt_sam, var.equal = var.equal, paired = FALSE)
  SG_Pval_FC <- aLL_Pval_FC[which(aLL_Pval_FC$Padj < 0.05 & aLL_Pval_FC$log2FC > FC), ]
  SG_Exp <- exp_mat[row.names(SG_Pval_FC), ]
  SG_list <- list(aLL_Pval_FC, SG_Pval_FC, SG_Exp)
  names(SG_list) <- c('AGs_pfc', 'SGs_pfc', 'SGs_exp')
  return(SG_list)
}
##############################################