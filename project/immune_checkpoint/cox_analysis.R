---
  ###  title: "R code for univariate cox model"
  ###  author: "JunShang"
  ###  date: "20161112"
  
library(survival)
library(BhGLM)
library(glmnet)

### t: time of metastasis
### d: status of metastasis
### x.c: clinical factors
### x.g: gene expression


# ******************************************************************************
# clinical factors
cnam <- "LUAD"
dir_GDC <- "/Users/stead/Desktop/GDC/"
gnam <- c("ENSG00000260695.1",  "ENSG00000248608.2",  "ENSG00000168878.15",  "ENSG00000203878.10", 
          "ENSG00000112175.7",  "ENSG00000166961.13")
goffm <- c("RP11-513N24.11", "RP11-206P5.2", "SFTPB", "CHIAP2", "BMP5", "MS4A15")
rt_cox <- read.table(file = paste(dir_GDC, cnam , "/cox/data_cox_t.tsv", sep = ""), header = T, as.is = T,  na.strings = "NA", row.names = NULL, sep = "\t")
rt_cox_M <- rt_cox[, c("row.names", "RFS_time", "S_stat", "Gender", "M", "N", "T", "pathologic_stage", 
                       "Age", "tobacco", gnam)]
colnames(rt_cox_M) <- c("sampleID", "time", "status", "gender", "M", "N", "T", "stage", 
                        "age", "smoke", goffm)
x.c <- rt_cox_M[which(rt_cox_M[, "status"] != "NA"), ]

#a <- c(sample(length(row.names(x.c))))
#data.T.test <- x.c[a[1: (length(x.c[, 1]) %/% 2)], ] ##
#dim(data.T.test)
#data.T.validate <- x.c[a[(length(x.c[, 1]) %/% 2 + 1) : length(x.c[, 1])], ]
#dim(data.T.validate)
#x.c <- data.T.validate

t <- x.c[, "time"]
d <- x.c[, "status"]
y <- Surv(t, d)


### docox a the function counted the result of cox
docox <- function(x){
  xc = x.c[, x]
  x1 = xc[, 3:5]
  x2 = xc[, 1:2]
  apply(x2, 2, table)
  ### can rescale you parameter, when you data have NA or NULL, you can use it
  xc = covariates(x1, x2, rescale = T, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  return(fc)
}

coxout <- function(x){
  k = c("gender",	"stage", "smoke", "age")
  NewV = append(k, goffm[x], after = 4)  ###add variable value behind the a. the after = 7 set which site of a begin to add variable value
  FC = docox(NewV)
  tt = summary(FC)
  rt_MVM <- tt[7:8] ### chose which results of summary(FC) will output
  tiff(file = paste(goffm[x], ".tiff", sep = ""), width = 600, height = 400)
  par(mar= c(4,15,8,8))
  plot.bh(FC, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE,
          cex = 2, cex.main = 2, cex.var= 2, cex.pts = 1, lwd= 3, OR = T, main = goffm[x])
  dev.off()
  return(rt_MVM)
}
coxout(5)


#b <- c(colnames(x.c[15: length(colnames(x.c))])) ### the variable value will be add to a 
count_cox_tab <- function(x){
  data_MVM <- coxout(x)
  MVM <- data.frame(matrix(unlist(data_MVM), nrow=5, byrow=F),stringsAsFactors=FALSE)
  row.names(MVM) <- c("smoke", "age", goffm[x], "gender",	"stage")
  colnames(MVM) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", " exp(coef)", 
                     "exp(-coef)", "lower0.95", "upper0.95")
  write.table(MVM, file = paste(goffm[x], "_MVM", ".tsv", sep = ""), row.names = T, col.names = T, sep = "\t")
}

x <- c(1:6)
lapply(x, count_cox_tab)

### u cox model for the clinical factors
### x is used to choose  gene, you can get the univariate analysis results of all parameters
x <- 6
uvm_count <- function(i, x){
  k = c("gender",	"stage", "smoke", "age")
  NewV = append(k, goffm[x], after = 4) 
  xc = x.c[, NewV]
  x1 = xc[, 3:5]
  x2 = xc[, 1:2]
  apply(x2, 2, table)
  ### can rescale you parameter, when you data have NA or NULL, you can use it
  xc = covariates(x1, x2, rescale = T, fill.missing = T)
  fc = coxph(y ~xc[, i], data = xc)
  tt <- summary(fc)
  rt_UVM <- tt[7:8]
  return(rt_UVM)
}

i <- c(1: 5)
data_UVM <- lapply(i, uvm_count, x)
UVM <- data.frame(matrix(unlist(data_UVM), nrow=5, byrow=T),stringsAsFactors=FALSE)
row.names(UVM) <- c("smoke", "age", goffm[x], "gender",	"stage")
colnames(UVM) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", " exp(coef)", 
                   "exp(-coef)", "lower0.95", "upper0.95")
