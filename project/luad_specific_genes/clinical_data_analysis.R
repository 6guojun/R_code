setwd("/Users/stead/Desktop/mir-370/clinical")
cli_group <- function(x){
  
  x$T[grep("Tis", x$T)] <- "1"
  x$T[grep("T0", x$T)] <- "1"
  x$T[grep("T1", x$T)] <- "1"
  x$T[grep("T2", x$T)] <- "1"
  x$T[grep("T3", x$T)] <- "2"
  x$T[grep("T4a", x$T)] <- "2"
  x$T[grep("T4b", x$T)] <- "2"
  x$T[grep("T4", x$T)] <- "2"
  x$T[grep("TX", x$T)] <- ""
  x$T[grep("Discrepancy", x$T)] <- ""
  
  
  x$N[grep("N0", x$N)] <- "1"
  x$N[grep("N1", x$N)] <- "2"
  x$N[grep("N2a", x$N)] <- "2"
  x$N[grep("N2b", x$N)] <- "2"
  x$N[grep("N2c", x$N)] <- "2"
  x$N[grep("N2", x$N)] <- "2"
  x$N[grep("N3", x$N)] <- "2"
  x$N[grep("NX", x$N)] <- ""
  x$N[grep("Discrepancy", x$N)] <- ""
  
  x$M[grep("M0", x$M)] <- "1"
  x$M[grep("M1", x$M)] <- "2"
  x$M[grep("MX", x$M)] <- ""
  x$M[grep("Discrepancy", x$M)] <- ""
  
  x$stage[grep("Stage IVA", x$stage)] <- "2"
  x$stage[grep("Stage IVB", x$stage)] <- "2"
  x$stage[grep("Stage IVC", x$stage)] <- "2"
  x$stage[grep("Stage III", x$stage)] <- "2"
  x$stage[grep("Stage II", x$stage)] <- "1"
  x$stage[grep("Stage I", x$pstage)] <- "1"
  x$stage[grep("[Discrepancy]", x$stage)] <- ""
  
  return(x)
}

#gnam <- "MIMAT0026483"
gnam <- "MIMAT0004680"
cnam <- "COAD"
miRNA_id_T <- read.table(file = "/Users/stead/Desktop/GDC/miRNA_id_transform.txt", header = T, sep = "\t")
rt <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnam, "/raw_data/", cnam, "_miRNA_HiSeq_gene", sep = ""), header = T, row.names = 1, sep = "\t")
colnames(rt) <- colnames(rt)[-1]
rt <- t(rt)
exp_rt <- cbind(row.names(rt), rt)
cli_rt <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnam , "/clinical_data/", cnam, "_clinicalMatrix", sep = ""), header = T, as.is = T,  na.strings = "NA", row.names = NULL, sep = "\t")
cli_rt[, "sampleID"] <- gsub("-", ".", cli_rt[, "sampleID"])
pre_data <- merge(cli_rt, exp_rt, by = 1)
sur_rt <- pre_data[c("pathologic_T", "pathologic_M", "pathologic_N", "pathologic_stage", gnam)]
### you can change the colnames or not
nam <-  as.character(miRNA_id_T[which(miRNA_id_T[, 1] == gnam), 2])
colnames(sur_rt) <- c("T", "M", "N", "stage", nam)

data_cli <- cli_group(sur_rt)
data_cli[is.na(data_cli)] <- ""
data_cli[, nam] <- as.numeric(as.character(data_cli[, nam]))
data_cli[, nam][is.na(data_cli[, nam])] <- 0


mak_cli <- function(data_cli, k){
  data_t <- data_cli[which(data_cli[, k] != ""), ]
  if(is.null(k)||k == ""){
    return(NULL)
  } else {
    data_e_t <- data_t[, c(k, nam)]
    data_e_t <- data_e_t[which(data_e_t[, k] != ""), ]
    return(data_e_t)
  }
}

count_cli <- function(data_cli, k){
  data_e_t <- mak_cli(data_cli, k)
  T_value <- t.test(data_e_t[, nam]~data_e_t[, k], paired = FALSE)
  tiff(file = paste(cnam, "_",  nam, "_", k, ".tiff", sep = ""), width = 512, height = 512)
  txt <- paste("p-value=", round(T_value$p.value[1], digits = 4), sep = "")
  par(mar=c(5, 5, 5, 2 + 0.1))
  plot(as.factor(data_e_t[, k]), data_e_t[, nam], ylab = "expression_level", main = txt, 
       col = c("blue", "red"), xlab = paste(cnam, k), cex.axis = 2, cex.main = 2, cex.lab = 2)
  dev.off()
}


count_cli(data_cli, "T")
count_cli(data_cli, "N")
count_cli(data_cli, "M")
count_cli(data_cli, "stage")




