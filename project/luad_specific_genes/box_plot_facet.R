library(ggplot2)

a <- c("ENSG00000019169.10", "ENSG00000122852.13", "ENSG00000168484.11", "ENSG00000185303.14", 
       "ENSG00000203878.10", "ENSG00000231322.4", "ENSG00000260695.1", "ENSG00000248608.2",
       "ENSG00000168878.15", "ENSG00000203878.10", "ENSG00000112175.7", "ENSG00000166961.13")

b <- c("MARCO", "SFTPA1", "SFTPC", "SFTPA2", "CHIAP2", "RPL13AP17", "RP11-513N24.11", "RP11-206P5.2", 
       "SFTPB", "CHIAP2", "BMP5", "MS4A15")


cnam <- "LUAD"
rt <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnam, "/raw_data/", cnam, "_mRNA_FPKM.txt", sep = ""), header = T, sep = "\t")

splite_data <- function(rt, k){
  ### splite the tumor and normal base on k which is a character 
  
  keepRow <- c(1: nrow(rt))            
  keepCol <- seq(1, ncol(rt)) ## product the data of read per millions
  data <- rt[keepRow, keepCol]
  data <- log2(data + 0.001)
  dimnames <- list(rownames(data), colnames(data))
  data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames) 
  
  G.R <- sapply(strsplit(colnames(data), "\\."), "[", 4 ) # "["(a, 4) almost a[4]
  
  
  if(is.null(k)||k == ""){ 
    ### you can choose the output data base on k
    return(NULL)
  } else {
    data.index <- grep(k, G.R)
    data_out <- data[, data.index]
    return(data_out) 
  }
}

mRNA_T_out <- splite_data(rt, "0..")
mRNA_N_out <- splite_data(rt, "1..")

mak_merge <- function(k, mRNA_T_out, mRNA_N_out){
  rt_gnam_T <- data.frame(t(rbind("Tumor", mRNA_T_out[a[k], ])), stringsAsFactors = F)
  rt_gnam_N <- data.frame(t(rbind("Normal", mRNA_N_out[a[k],])), stringsAsFactors = F)
  rt_gnam <- cbind(b[k], rbind(rt_gnam_N, rt_gnam_T))
  return(rt_gnam)
}

t <- c(1:length(a))
data_M <- lapply(t, mak_merge, mRNA_T_out, mRNA_N_out)

 
