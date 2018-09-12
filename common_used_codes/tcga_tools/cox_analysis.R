library(survival)
library(BhGLM)
library(glmnet)

uvm_count <- function(unam, rt_uvm){
  #rt_uvm contain OS_Time, OS_Status and other element 
  t <- rt_uvm[, "OS_Time"]
  d <- rt_uvm[, "OS_Status"]
  y <- Surv(t, d)
  
  x =  as.numeric(as.character(rt_uvm[, unam]))
  fc = coxph(y ~ x)
  tt <- summary(fc)
  uvm_out <- data.frame(cbind(tt$coefficients, tt$conf.int), stringsAsFactors = FALSE)[, c("coef", "exp.coef.", "lower..95", "upper..95", "Pr...z.." )]
  colnames(uvm_out) <- c("coef", "exp_coef", "lower_95%CI", "upper_95%CI", "pvalue")
  uvm_out <- data.frame(cbind(row.names(uvm_out), uvm_out))
  colnames(mvm_out)[1] <- 'variable'
  return(uvm_out)
}


#uvm_list <- lapply(unam, uvm_count, rt_uvm)
#data_uvm <- do.call(rbind, uvm_list)


mvm_count <- function(con_nam, log_nam, tnam, rt_esc){
  t <- rt_esc[, "OS_Time"]
  d <- rt_esc[, "OS_Status"]
  y <- Surv(t, d)
  
  CountCoxph <- function(on_nam, log_nam, tnam, rt_esc){
    if(is.na(log_nam)){
      xc_con = rt_esc[, con_nam]
      ### can rescale you parameter, when you data have NA or NULL, you can use it
      fc = coxph(y ~ ., data = xc_con) 
      return(fc)
    } else if (is.na(con_nam)){
      return(NULL)
    } else {
      xc_con = rt_esc[, con_nam]
      xc_log = rt_esc[, log_nam]
      ### can rescale you parameter, when you data have NA or NULL, you can use it
      xc = covariates(xc_con, xc_log, con.rescale = T, fill.missing = T)
      fc = coxph(y ~ ., data = xc)
      return(fc)
    }
  }
  
  fc <- CountCoxph(on_nam, log_nam, tnam, rt_esc)
  tt <- summary(fc)
  mvm_out <- data.frame(cbind(tt$coefficients, tt$conf.int), stringsAsFactors = FALSE)[, c("coef", "exp.coef.", "lower..95", "upper..95", "Pr...z.." )]
  colnames(mvm_out) <- c("coef", "HR", "lower_95%CI", "upper_95%CI", "pvalue")
  mvm_out <-data.frame(cbind(row.names(mvm_out), mvm_out))
  colnames(mvm_out)[1] <- 'variable'  

  pdf(file = paste(tnam, "_cox.pdf",  sep = ""), width = 10, height = 10)
  par(cex = 3)
  plot.bh(fc, threshold = 0.05, show.all.vars = T, col.pts = c("red", "black"), gap = 0, show.pvalues = TRUE,
          cex.var= 2, cex.pts = 0.5, OR = T, lwd = 4)
  dev.off()
  return(mvm_out)
}

#uvm_genes_pos <- c(7296, 24373, 16987)
#rt_mvm_d <-rt_T_m[uvm_genes_pos, ]
#mnam <- rt_exp$Ensembl_ID[uvm_genes_pos] 

#rt_mvm <- data.frame(cbind(t(rt_mvm_d), rt_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)
#colnames(rt_mvm) <- c(mnam, "OS_Status", "OS_Time")

