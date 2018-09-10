###usage: CountTimeROC(rt_mvm_sur_min, 'OS_Time', 'OS_Status', 'risk_score', predict.time = 365, lambda = 0.05, 1)

library(survival)
library(survivalROC)

CountTimeROC <- function(rt_roc, OS_Time, OS_Status, OS_marker){
  risk_roc1 = survivalROC(Stime = rt_roc[, OS_Time], status = rt_roc[, OS_Status], marker = rt_roc[, OS_marker], predict.time = 365, lambda = 0.05)
  risk_roc2 = survivalROC(Stime = rt_roc[, OS_Time], status = rt_roc[, OS_Status], marker = rt_roc[, OS_marker], predict.time = 1095, lambda = 0.05)
  risk_roc3 = survivalROC(Stime = rt_roc[, OS_Time], status = rt_roc[, OS_Status], marker = rt_roc[, OS_marker], predict.time = 1825, lambda = 0.05)
  
  pdf(file = paste(OS_marker, '_time_roc.pdf', sep = ""), 10, 10)
  par(cex = 2)
  plot(risk_roc1$FP, risk_roc1$TP, type = "l", xlim=c(0, 1), ylim=c(0, 1), col = 'red', lwd = 3, pch = 19, 
       xlab=paste( "False positive rate"),
       ylab= "True positive rate", main ="Method = NNE")
  lines(risk_roc2$FP, risk_roc2$TP, col = "blue", lwd = 3)
  lines(risk_roc3$FP, risk_roc3$TP, col = "cyan3", lwd = 3)
  
  legend("bottomright", legend=c(paste("Year = 1", paste(' AUC = ', round(risk_roc1$AUC, 3), sep = "")), 
                                 paste("Year = 3", paste(' AUC = ', round(risk_roc2$AUC, 3), sep = "")), 
                                 paste("Year = 5", paste(' AUC = ', round(risk_roc3$AUC, 3), sep = ""))),
         col=c('red', 'blue', 'cyan3'), lty = 1, lwd = 3,  cex = 1)
  
  abline(0, 1)
  dev.off()
}
