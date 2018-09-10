#################################################
###two table you should prepare. The colname is samples' name and the rowname is the gene name in expression table.
###group.txt file contain two column which contain samples' name in accord with expression table. 
###the other column show divide the samples into different group 
###you can use this function to get FC and pValue of the different expression genes, 
###1. DEGPvalFC <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "DEGPvalFC")
###2. DEGExp <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "DEGExp")
###3. AllPvalFC <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "AllPvalFC")
###4. DEGTest("H9D0", "H9D2", rt, rt_sam, FC,  var.equal = FALSE, paired = FALSE, meas = "")
###20170907
###JunShang
###Email:shangjunv@163.com
##################################################



##################################################
##prepared the table of gene expression for test##
#setwd("/Users/stead/Documents/SourceTree/R/RNA_seq/DEG_analysis")
#rt <- read.table(file ="WJC_H9_FPKM.txt")
#rt <- rt[apply(rt, 1, mean) > 0, ]
#rt_sam <- read.table(file = "group.txt", header = TRUE, sep = "\t", row.name = NULL, stringsAsFactors = FALSE)
##################################################

##################################################
###### main function #############################

print("your data must be concert to log format")
DEGTest <- function(GP1, GP2, rt, rt_sam, FC, var.equal = var.equal, paired = FALSE, meas = c("DEGPvalFC", "DEGExp", "AllPvalFC")){
  #you can get one table which contain information of pvalue and fold change
  #you can aslo get one expression table of DEG 
  #you can choice achive two table as above
  
  DrawMat <- function(GP1, GP2, rt, rt_sam){
    #GP2 and GP2 are the name of your group
    #rt is the expression table 
    #rt_sam contain samples name and its group name
    #this function will get a new expression table which show the group information in colnames
    
    M_Position <- match(colnames(rt), rt_sam$samples, nomatch = NA)
    Group_Nam <- rt_sam$group[M_Position]
    colnames(rt) <- paste(colnames(rt), Group_Nam, sep = "_")
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
    if(all(x[1] == (c(x, y)))){
      return(NA)
    } else{
      Obj <- t.test(x, y, var.equal = var.equal, paired = FALSE)
      Pval <- Obj$p.value
      Padj <-  p.adjust(Pval, method = "fdr")
      Cval <- c(Pval, Padj)
      return(Cval) 
    }
  }
  
  GetFC <- function(Gnam, GP1, GP2, data){
    #get the fold changeo of one gene between two group
    
    data_Gnam <- data[Gnam, ]
    GPx <- c(grep(GP1, colnames(data_Gnam)))
    GPy <- c(grep(GP2, colnames(data_Gnam)))
    x <- as.numeric(data_Gnam[, GPx])
    y <- as.numeric(data_Gnam[, GPy])
    log2FC <-  mean(y) - mean(x)
    return(log2FC)
  }
  
  MainTest <- function(GP1, GP2, rt, rt_sam, var.equal = var.equal, paired = FALSE){
    #get one table which contain information of pvalue and fold change
    
    rt_data <- DrawMat(GP1, GP2, rt, rt_sam)
    Gnam <- row.names(rt_data)
    Cval_L <- lapply(Gnam, GetPval, GP1, GP2, rt_data, var.equal = var.equal, paired = FALSE)
    data_Pval <- do.call(rbind, Cval_L)
    Log2FC_L <- lapply(Gnam, GetFC, GP1, GP2, rt_data)
    data_FC <- do.call(rbind, Log2FC_L)
    data_Pval_FC <- data.frame(cbind(data_Pval, data_FC), stringsAsFactors = FALSE)
    row.names(data_Pval_FC) <- row.names(rt_data)
    colnames(data_Pval_FC) <- c("Pval", "Padj", "log2FC")
    return(data_Pval_FC)
  }
  
  ###get the Pval and FC table
  ALL_Pval_FC <- MainTest(GP1, GP2, rt, rt_sam, var.equal = var.equal, paired = FALSE)
  DEG_Pval_FC <- ALL_Pval_FC[which(ALL_Pval_FC$Padj < 0.05 & abs(ALL_Pval_FC$log2FC) > FC), ]
  data_G1_G2 <- DrawMat(GP1, GP2, rt, rt_sam)
  DEG_Table <- data_G1_G2[row.names(DEG_Pval_FC), ]
  if (meas == ""|is.null(meas) == TRUE|all(meas == c("DEGPvalFC", "DEGExp", "AllPvalFC")) == TRUE){
    write.table(DEG_Pval_FC, file = paste(GP1, "_", GP2, "_DEG_", "Pval", "_", "FC", ".txt", sep = ""), sep = "\t")
    write.table(DEG_Table, file = paste(GP1, "_", GP2, "_DEG_", "_Exp_", ".txt", sep = ""), sep = "\t") 
    write.table(ALL_Pval_FC, file = paste(GP1, "_", GP2, "_ALL", "_Pval_", "FC", ".txt", sep = ""), sep = "\t")
  } else if (meas == "DEGPvalFC"){
    return(DEG_Pval_FC)
  } else if (meas == "DEGExp"){
    return(DEG_Table)
  } else if (meas == "AllPvalFC") {
    return(ALL_Pval_FC)
  } else {
    print("if you write parameter, meas must be PvalFC or DEGExp")
  }
}
##############################################
