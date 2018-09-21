CountCindex <- function(rt_exp_sur, Gnam){
  fit <- coxph(Surv(OS_Time, OS_Status) ~ rt_exp_sur[, Gnam], data = rt_exp_sur)
  cindex <- concordance.index(predict(fit), surv.time = rt_exp_sur$OS_Time, surv.event = rt_exp_sur$OS_Status, method = "noether")
  c_index <- round(cindex$c.index, digits = 4)
  pvalue <- round(cindex$p.value, digits = 4)
  print(paste('c-index: ', c_index, '/', 'pvalue: ', pvalue, sep = ''))
  return(c_index)
}

