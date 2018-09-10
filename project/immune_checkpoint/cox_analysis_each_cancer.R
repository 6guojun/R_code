###  title: "R code for univariate cox model"
###  author: "JunShang"
###  date: "20161112"

library(survival)
library(BhGLM)
library(clusterProfiler)

### t: time of metastasis
### d: status of metastasis
### x.c: clinical factors
### x.g: gene expression
setwd("/Users/stead/Desktop/PD-L1_and_TMI_type/cox_analysis")

# ******************************************************************************
# clinical factors
gnam <- "ENSG00000188389"
###expression data
source("/Users/stead/Desktop/PD-L1_and_TMI_type/scripts/get_exp_sur_cli_merge_data.R")

#############################################LUAD####################################################
###get merge data

LUADCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic", "number_pack_years_smoked", "race.demographic",
                             "pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "smoke", "race", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$race <- factor(rt_esc$race, labels=c(1:5))
  rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 4, 5))
  rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(5, 1, 1, 1, 1, 1, 1, 2, 2, 2))
  rt_esc[, c("gender", "race", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "race", "T", "N", "M", "stage")], 
                                                                           ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "smoke", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 5)), ]
}


LUADCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4, 6)]
  x2 = rt_esc_d[, c(5, 7)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf",  sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
}

LUADMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = LUADCliSort(gnam, cancer, data_list)
  LUADCox(rt_esc_d, cancer)
}

LUAD_list <- get_exp_sur_cli("LUAD")
cox_table_list <- lapply(gnam, LUADMainFun, "LUAD", LUAD_list)
cox_LUAD <- do.call(rbind, cox_table_list)
write.table(cox_LUAD, file = paste("LUAD", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################UVM####################################################

UVMCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses",  "new_tumor_event_after_initial_treatment")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage", "recurrence")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
 
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3))
  rt_esc$N <- factor(rt_esc$N, labels=c(2, 1, 2))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2))
  rt_esc$recurrence <- factor(rt_esc$recurrence, labels=c(3, 1, 2))
  rt_esc[, c("gender", "T", "N", "M", "stage", "recurrence")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage", "recurrence")], 
                                                                           ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-c(which(rt_esc$stage == 3), which(rt_esc$recurrence == 3)), ]
}

UVMCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
   
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf",  sep = ""), width = 10, height = 10)
  par(cex = 3, cex.axis = 5)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

UVMMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = UVMCliSort(gnam, cancer, data_list)
  UVMCox(rt_esc_d, cancer)
}

UVM_list <- get_exp_sur_cli("UVM")
cox_table_list <- lapply(gnam, UVMMainFun, "UVM", UVM_list)
cox_UVM <- do.call(rbind, cox_table_list)
write.table(cox_UVM, file = paste("UVM", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################BRCA####################################################

BRCACliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5))
  rt_esc$N <- factor(rt_esc$N, labels=c(rep(1, 4), rep(2, 11), 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(1, 2, 3, 1))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 3)), ]
}

BRCACox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

BRCAMainFun <- function(gnam, cancer, rt_T_m){
  rt_esc_d = BRCACliSort(gnam, cancer, rt_T_m)
  BRCACox(rt_esc_d, cancer)
}

BRCA_list <- get_exp_sur_cli("BRCA")
cox_table_list <- lapply(gnam, BRCAMainFun, "BRCA", BRCA_list)
cox_BRCA <- do.call(rbind, cox_table_list)
write.table(cox_BRCA, file = paste("BRCA", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################ESCA####################################################

ESCACliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(6, 1, 2, 3, 4, 5, 5))
  rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 3)), ]
}

ESCACox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

ESCAMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = ESCACliSort(gnam, cancer, data_list)
  ESCACox(rt_esc_d, cancer)
}

ESCA_list <- get_exp_sur_cli("ESCA")
cox_table_list <- lapply(gnam, ESCAMainFun, "ESCA", ESCA_list)
cox_ESCA <- do.call(rbind, cox_table_list)
write.table(cox_ESCA, file = paste("ESCA", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################HNSC####################################################

HNSCCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(6, 1, 2, 3, 4, 5, 5, 5, 6))
  rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 2, 2, 2, 2, 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 1, 1, 2, 2, 2, 2))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 3)), ]
}

HNSCCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

HNSCMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = HNSCCliSort(gnam, cancer, data_list)
  HNSCCox(rt_esc_d, cancer)
}

HNSC_list <- get_exp_sur_cli("HNSC")
cox_table_list <- lapply(gnam, HNSCMainFun, "HNSC", HNSC_list)
cox_HNSC <- do.call(rbind, cox_table_list)
write.table(cox_HNSC, file = paste("HNSC", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)


#############################################KIRP####################################################

KIRPCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5))
  rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels = c(3, 1, 1, 2, 2))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 3)), ]
}

KIRPCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

KIRPMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = KIRPCliSort(gnam, cancer, data_list)
  KIRPCox(rt_esc_d, cancer)
}

KIRP_list <- get_exp_sur_cli("KIRP")
cox_table_list <- lapply(gnam, KIRPMainFun, "KIRP", KIRP_list)
cox_KIRP <- do.call(rbind, cox_table_list)
write.table(cox_KIRP, file = paste("KIRP", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################KIRC####################################################

KIRCCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4))
  rt_esc$N <- factor(rt_esc$N, labels=c(1, 2, 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels = c(3, 1, 1, 2, 2))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 3)), ]
}

KIRCCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

KIRCMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = KIRCCliSort(gnam, cancer, data_list)
  KIRCCox(rt_esc_d, cancer)
}

KIRC_list <- get_exp_sur_cli("KIRC")
cox_table_list <- lapply(gnam, KIRCMainFun, "KIRC", KIRC_list)
cox_KIRC <- do.call(rbind, cox_table_list)
write.table(cox_KIRC, file = paste("KIRC", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################LAML####################################################

LAMLCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5))
  rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, 2, 2, 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 3))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(5, 1, 2, 3, 4))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 5)), ]
}

LAMLCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

LAMLMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = LAMLCliSort(gnam, cancer, data_list)
  LAMLCox(rt_esc_d, cancer)
}

LAML_list <- get_exp_sur_cli("LAML")
cox_table_list <- lapply(gnam, LAMLMainFun, "LAML", LAML_list)
cox_LAML <- do.call(rbind, cox_table_list)
write.table(cox_LAML, file = paste("LAML", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)


#############################################LGG####################################################

LGGCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","neoplasm_histologic_grade")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "grade")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$grade <- factor(rt_esc$grade, labels=c(3, 1, 2))
  rt_esc[, c("gender", "grade")] <- apply(as.matrix(rt_esc[, c("gender", "grade")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "grade")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$grade == 3)), ]
}

LGGCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

LGGMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = LGGCliSort(gnam, cancer, data_list)
  LGGCox(rt_esc_d, cancer)
}

LGG_list <- get_exp_sur_cli("LGG")
cox_table_list <- lapply(gnam, LGGMainFun, "LGG", LGG_list)
cox_LGG <- do.call(rbind, cox_table_list)
write.table(cox_LGG, file = paste("LGG", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################OV####################################################

OVCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "clinical_stage")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "clinical_stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$clinical_stage <- factor(rt_esc$clinical_stage, labels=c(3, 1, 1, 1, 1, 2, 2, 2, 2))
  rt_esc[, "clinical_stage"] <- apply(as.matrix(rt_esc[, "clinical_stage"], 
                                                    ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "clinical_stage")]
  rt_esc_d <- rt_esc_d[-which(rt_esc$clinical_stage == 3 ), ]
}

OVCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, 5]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

OVMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = OVCliSort(gnam, cancer, data_list)
  OVCox(rt_esc_d, cancer)
}


OV_list <- get_exp_sur_cli("OV")
cox_table_list <- lapply(gnam, OVMainFun, "OV", OV_list)
cox_OV <- do.call(rbind, cox_table_list)
write.table(cox_OV, file = paste("OV", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)

#############################################SKCM####################################################

SKCMCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","pathologic_T", "pathologic_N", "pathologic_M", "tumor_stage.diagnoses")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "T", "N", "M", "stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$gender <- factor(rt_esc$gender, labels=c(1, 2))
  rt_esc$T <- factor(rt_esc$T, labels=c(5, 1, 2, 2, 2, 3, 3, 4, 4, 4))
  rt_esc$N <- factor(rt_esc$N, labels=c(3, 1, rep(2, 8), 3))
  rt_esc$M <- factor(rt_esc$M, labels=c(3, 1, 2, 2, 2))
  rt_esc$stage <- factor(rt_esc$stage, labels=c(3, 3, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2))
  rt_esc[, c("gender", "T", "N", "M", "stage")] <- apply(as.matrix(rt_esc[, c("gender", "T", "N", "M", "stage")], 
                                                                   ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "gender", "stage")]
  rt_esc_d <- rt_esc_d[-(which(rt_esc$stage == 3)), ]
}

SKCMCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5, 6)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

SKCMMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = SKCMCliSort(gnam, cancer, data_list)
  SKCMCox(rt_esc_d, cancer)
}

SKCM_list <- get_exp_sur_cli("SKCM")
cox_table_list <- lapply(gnam, SKCMMainFun, "SKCM", SKCM_list)
cox_SKCM <- do.call(rbind, cox_table_list)
write.table(cox_SKCM, file = paste("SKCM", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)


#############################################UCEC####################################################

UCECCliSort <- function(gnam, cancer, data_list){
  rt_T_m <- data_list[[1]]
  rt_S_m <- data_list[[2]]
  rt_C_m <- data_list[[3]]
  
  rt_esc <- cbind(t(rt_T_m[grep(gnam, row.names(rt_T_m)), ]), rt_S_m[, c("X_OS_IND", "X_OS")], 
                  rt_C_m[, c("age_at_initial_pathologic_diagnosis", "gender.demographic","clinical_stage")])
  colnames(rt_esc) <- c(gnam, "status", "time",  "age", "gender", "clinical_stage")
  write.table(rt_esc, file = paste(cancer, "_exp_cli.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  ###convert clinical element to numeric
  rt_esc$clinical_stage <- factor(rt_esc$clinical_stage, labels=c(rep(1, 4), rep(1, 3), rep(2, 6), rep(2, 3)))
  rt_esc[, "clinical_stage"] <- apply(as.matrix(rt_esc[,  "clinical_stage"], ncol = 6, byrow = TRUE), 2, function(x){as.numeric(x)})
  rt_esc_d <- rt_esc[, c(gnam, "status", "time",  "age", "clinical_stage")]
}

UCECCox <- function(rt_esc_d, cancer){
  t <- rt_esc_d[, "time"]
  d <- rt_esc_d[, "status"]
  y <- Surv(t, d)
  
  x1 = rt_esc_d[, c(1, 4)]
  x2 = rt_esc_d[, c(5)]
  apply(x2, 2, table)
  xc = covariates(x1, x2, con.rescale = TRUE, cat.center = TRUE, fill.missing = T)
  fc = coxph(y ~ ., data = xc)
  tt <- summary(fc)
  pdf(file = paste(cancer,  "_", gnam, "_cox.pdf", sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE, 
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  cox_table <- tt$coefficients
}

UCECMainFun <- function(gnam, cancer, data_list){
  rt_esc_d = UCECCliSort(gnam, cancer, data_list)
  UCECCox(rt_esc_d, cancer)
}

UCEC_list <- get_exp_sur_cli("UCEC")
cox_table_list <- lapply(gnam, UCECMainFun, "UCEC", UCEC_list)
cox_UCEC <- do.call(rbind, cox_table_list)
write.table(cox_UCEC, file = paste("UCEC", "_cox.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE)
