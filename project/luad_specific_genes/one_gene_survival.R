library(ggplot2)
library(survival)
library(BhGLM)
library(glmnet)
library(gplots)
library(reshape)
### t: time of metastasis
### d: status of metastasis
### x.c: clinical factors
### x.g: gene expression


# ******************************************************************************
# clinical factors
cnam <- "LUAD"
dir_GDC <- "/Users/stead/Desktop/GDC/"
#gnam <- c("ENSG00000248608.2",  "ENSG00000168878.15", "ENSG00000112175.7")
#goffm <- c("RP11-206P5.2", "SFTPB", "BMP5")
gnam <- c("ENSG00000260695.1",  "ENSG00000248608.2",  "ENSG00000168878.15",  "ENSG00000203878.10", 
          "ENSG00000112175.7",  "ENSG00000166961.13")
goffm <- c("RP11-513N24.11", "RP11-206P5.2", "SFTPB", "CHIAP2", "BMP5", "MS4A15")

rt_cox <- read.table(file = paste(dir_GDC, cnam , "/cox/data_cox_t.tsv", sep = ""), header = T, as.is = T,  na.strings = "NA", row.names = NULL, sep = "\t")
rt_cox_M <- rt_cox[, c("row.names", "S_time", "S_stat", "RFS_time", "RFS_stat", "Gender", "M", "N", "T", "pathologic_stage", 
                       "Age", "tobacco", gnam)]
colnames(rt_cox_M) <- c("sampleID", "time", "status", "RFS_time", "RFS_stat", "gender", "M", "N", "T", "stage", 
                        "age", "smoke", goffm)
x.c <- rt_cox_M[which(rt_cox_M[, "status"] != "NA"), ]

count_sur <- function(y, x){
  a <- x[, y] < median(x[, y])
  diff <- survdiff(Surv(time, status) ~a, data = x)
  pValue <- 1-pchisq(diff$chisq,df=1)
  pValue <- round(pValue, 5)
  fit <- survfit(Surv(time, status) ~ a, data = x)
  summary(fit)   
  filenames = y
  pdf(file = paste(filenames,  ".pdf", sep = ""))
  par(mar = c(5, 5, 5, 3))
  plot(fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
       main = filenames)
  legend(0, 0.17, c("highrisk", "lowrisk"), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
  legend("topright", paste("p-value:", pValue, sep=" "), cex = 1.5)
  dev.off() 
}
lapply(goffm, count_sur, x.c)
