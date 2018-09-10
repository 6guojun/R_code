### you shou prepare miRNA expression data, miRNA id transform data and clinical data
### this code help you to do boxplot betwwen Tumor and Normal
### 

library(survival)
cnam <- "COAD"
gnam <- "hsa-miR-503-5p"

miRNA_exp <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnam, "/raw_data/", cnam, "_miRNA_HiSeq_gene", sep = ""), header = T, row.names = NULL, sep = "\t")
colnames(miRNA_exp) <- colnames(miRNA_exp)[-1]
miRNA_id_T <- read.table(file = "/Users/stead/Desktop/GDC/miRNA_id_transform.txt", header = T, sep = "\t")

splite_mi_data <- function(miRNA_exp, miRNA_id_T, k){
  G.R <- sapply(strsplit(colnames(miRNA_exp), "\\."), "[", 4 ) # "["(a, 4) almost a[4]
    if(is.null(k)||k == ""){ 
      ### TCGA-xx-xxxx-xx-xxxx. The sample id is the condition choosing the data_out
      return(NULL)
    } else {
      data.index <- grep(k, G.R) ## get the tumor group index
      miRNA_exp_T <- miRNA_exp[, data.index]
      miRNA_exp_T <- cbind(miRNA_exp[, 1], miRNA_exp_T)
      
      miRNA_id_T <- merge(miRNA_id_T, miRNA_exp_T, by = 1)
      miRNA_nam <- miRNA_id_T[, c(-1, -(3:6))]
      row.names(miRNA_nam) <- miRNA_nam[, 1]
      miRNA_T_out <- data.frame(apply(miRNA_nam[, -1], 2, as.numeric), stringsAsFactors = F)
      row.names(miRNA_T_out) <- miRNA_nam[, 1]
      return(miRNA_T_out)
    }
  }

  

miRNA_T_out <- splite_mi_data(miRNA_exp, miRNA_id_T, "0.")
miRNA_N_out <- splite_mi_data(miRNA_exp, miRNA_id_T, "1.")

mak_box_data <- function(miRNA_T_out, miRNA_N_out, gnam){
  ### prapare the data unit for boxplot
  g_tumor <- rbind("Tumor", miRNA_T_out[which(rownames(miRNA_T_out) == gnam), ])
  g_tumor <- data.frame(t(g_tumor), stringsAsFactors = F)
  g_normal <- rbind("Normal", miRNA_N_out[which(rownames(miRNA_N_out) == gnam), ])
  g_normal <- data.frame(t(g_normal), stringsAsFactors = F)
  data_unit <- rbind(g_tumor, g_normal)
  colnames(data_unit) <- c("group", gnam)
  data_unit[, 2] <- as.numeric(data_unit[, 2])
  return(data_unit)
}
data_unit <- mak_box_data(miRNA_T_out, miRNA_N_out, gnam)

mak_boxplot <- function(data_unit, gnam){
  T_Value <- t.test(data_unit[, 2]~group, data = data_unit, paired = FALSE)
  tiff(file = paste(cnam, gnam, ".tiff", sep = "" ), width = 216, height = 512)
  txt <- paste("p-value=", round(T_Value$p.value[1], digits = 4), sep = "")
  p = ggplot(data = data_unit, aes(data_unit[, 1], data_unit[, 2])) +
    xlab(cnam) +  ylab("gene expression") +
    geom_boxplot(fill = c("blue", "red")) + ggtitle(txt) +
    theme_bw()+
    theme(plot.title = element_text(size=20, face="bold", hjust = 0.5), 
          axis.title.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          axis.text.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
          axis.title.y = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
          axis.text.y = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
          panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5))
  print(p)
  dev.off()   
}

mak_boxplot(data_unit, gnam)

colnames(miRNA_T_out) <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*", "\\1\\.\\2\\.\\3\\.\\4", colnames(miRNA_T_out))
colnames(miRNA_T_out) <- substring(colnames(miRNA_T_out), 1, 15)

exp_rt <- t(miRNA_T_out) 
exp_rt <- cbind(row.names(exp_rt), data.frame(exp_rt, stringsAsFactors = F))
cli_rt <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnam, "/clinical_data/", cnam,  "_clinicalMatrix", sep = "") , header = T, row.names = NULL, sep = "\t")
cli_rt[, "sampleID"] <- gsub("-", ".", cli_rt[, "sampleID"])
pre_data <- merge(cli_rt, exp_rt, by = 1)

sur_data <- function(gnam, data){
  gnam <- gsub("-", ".", gnam)
  sur_rt <- pre_data[c("X_EVENT", "X_OS", gnam)]
  sur_rt <- sur_rt[which(sur_rt$"X_EVENT" != "NA"), ]
  sur_rt <- sur_rt[which(sur_rt[, gnam] != "NA"), ]
  colnames(sur_rt) <- c("status", "time", gnam)
  return(sur_rt)
}

count_sur <- function(gnam, data){
  gnam <- gsub("-", ".", gnam)
  sur_rt <- sur_data(gnam, pre_data)
  b  <- median(sur_rt[, gnam])
  a <- sur_rt[, gnam] < median(sur_rt[, gnam])
  if(b == 0){
    return(NULL)
  }else {
    diff <- survdiff(Surv(time, status) ~a, data = sur_rt)
    pValue <- 1-pchisq(diff$chisq,df=1)
    pValue <- round(pValue, 5)
    fit <- survfit(Surv(time, status) ~ a, data = sur_rt)
    summary(fit)   
    pdf(file = paste(cnam, gnam, ".pdf", sep = ""))
    par(mar = c(5, 5, 5, 3))
    plot(fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
         main = paste(gnam, "(", cnam, ")", sep="")) 
    legend(0, 0.17, c("highExp", "lowExp"), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
    legend("topright", paste("p-value:", pValue, sep=" "), cex = 1.5)
    dev.off() 
  }
}
#gnam <- colnames(exp_rt)[-1]
count_sur(gnam, pre_data)



