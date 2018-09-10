### Filter the cancer tissue special genes. 
### we first count the FC >1, Mean1 - mean2
### you should prepare all cancers RNA-seq data (FPKM)
### 20170223 sj

library("ggplot2")
setwd("/Users/stead/Desktop/GDC/LUAD/ROC/tissue_special_gene")
k <- "0.."

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



merge_mean <- function(cnam, k){
  rt <- read.table(file = paste(file = "/Users/stead/Desktop/GDC/", cnam, "/raw_data/", cnam, "_mRNA_FPKM.txt", sep = ""), header = T, na.strings = NA, row.names = 1, sep = "\t")
  if(is.null(k)||k == ""){ 
    return(NULL)
  } else {
    data_out <- splite_data(rt, k)
    data_out_M <- data.frame(apply(data_out, 1, mean), stringsAsFactors = F) 
    data_out_M <- cbind(row.names(data_out_M), data_out_M)
    return(data_out_M)
  }
}


cnam <- c("LUAD", "BRCA", "BLCA", "ACC", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "LUSC", "OV", "PAAD", "PRAD", "PCPG", "READ", "SARC", "SKCM",
"STAD", "TGCT", "THCA", "THYM", "UCEC", "LAML", "LGG")
#b <- c("LUAD", "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUSC", "STAD", "THCA", "UCEC")

rt_T_M <- lapply(cnam, merge_mean, "0..")
#rt_N_M <-  lapply(b, merge_mean, "1..")

rt_T_M_M <- do.call("cbind", rt_T_M)
#rt_N_M_M <- do.call("cbind", rt_N_M)

rt_M <- rt_T_M_M[, which(colnames(rt_T_M_M) == "apply.data_out..1..mean.")]
rt_M_name <- rt_T_M_M[, which(colnames(rt_T_M_M) == "row.names(data_out_M)")]

#rt_M <- rt_N_M_M[, which(colnames(rt_N_M_M) == "apply.data_out..1..mean.")]
#rt_M_name <- rt_N_M_M[, which(colnames(rt_N_M_M) == "row.names(data_out_M)")]

colnames(rt_M) <- cnam
#colnames(rt_M) <- b
#write.table(rt_M, file = "normal_mean_data.txt", row.names = T, col.names = T, sep = "\t")
#rite.table(rt_M_name, file = "normal_mean_data_name.txt", row.names = T, col.names = T, sep = "\t")

rt_M_t <- data.frame(t(rt_M), stringsAsFactors = F)
#colnames(rt_M_t) <- c(rt_M_t[2, ])
#data <- data.frame(apply(rt_M_t[-c(1,2), ], 2, as.numeric), stringsAsFactors = F)
rt_M_C <- cbind(row.names(rt_M_t), rt_M_t)
colnames(rt_M_C)[1] <- "id"


draw_box <- function(data_out){
  ### id is the group of cancer and data_out[2] contain all expression of these cancer
  filename <- colnames(data_out)[2]
  pdf(file = paste(filename, ".pdf", sep = ""))
  p = qplot(id, data_out[,2], data = data_out, geom = "boxplot") + 
      theme(plot.title = element_text(size=100, face="bold", hjust = 0.5), 
            axis.title.x = element_text(size = 7, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
            axis.text.x = element_text(size = 7, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450))
  print(p)
  dev.off()   
} 

max_count <- function(x){
  filename <- colnames(rt_M_C)[x]
  t <- which(row.names(rt_M_C) == "LUAD")
  r <- which(row.names(rt_M_C) != "LUAD")
  k <- rt_M_C[r, x]
  
  if(as.character(rt_M_C[which.max(rt_M_C[, x]), 1]) == "LUAD" &&
    #all(rt_M_C[t, x]/(k + 0.0001) >= 2) == TRUE) 
    all((rt_M_C[t, x] - k >= 1) == TRUE)){
    data_out <- rt_M_C[, c(1, x)]
    write.table(data_out, file = paste(filename, ".tsv", sep = ""), col.names = T, row.names = T, sep = "\t")
    draw_box(data_out)
  }
} 



b <- c(2 : length(colnames(rt_M_C)))
lapply(b, max_count) 

