library(forestplot)

setwd('/Users/stead/Desktop/subtype_analysis/signature/mut_platform/LUAD/cox_results')
data_dir <- '/Users/stead/Desktop/subtype_analysis/signature/mut_platform/LUAD/'

###load data
load(paste(data_dir, 'TCGA', '/cox/', 'TCGA', '_uvm_cox.Rdata', sep = ''))
TCGA_uvm <- rbind(NA, data_uvm_cox)
TCGA_uvm$variable <- c('TCGA', 'risk_score', 'age', 'T2', 'T3', 'T4', 'N1-3', 'M1', 'stageII', 'stageIII', 'stageIV', 'recurrence', 'gender')

load(paste(data_dir, 'GSE30219', '/cox/', 'GSE30219', '_uvm_cox.Rdata', sep = ''))
GSE30219_uvm <- rbind(NA, data_uvm_cox)
GSE30219_uvm$variable <- c('GSE30219', 'risk_score', 'age' , 'T2', 'N1', 'recurrence', 'gender')

load(paste(data_dir, 'GSE31210', '/cox/', 'GSE31210', '_uvm_cox.Rdata', sep = ''))
GSE31210_uvm <- rbind(NA, data_uvm_cox)
GSE31210_uvm$variable <- c('GSE31210', 'risk_score',  'age', 'stageII', 'recurrence', 'gender')

load(paste(data_dir, 'GSE3141', '/cox/', 'GSE3141', '_uvm_cox.Rdata', sep = ''))
GSE3141_uvm <- rbind(NA, data_uvm_cox)
GSE3141_uvm$variable <- c('GSE3141', 'risk_score')

load(paste(data_dir, 'GSE81089', '/cox/', 'GSE81089', '_uvm_cox.Rdata', sep = ''))
GSE81089_uvm <- rbind(NA, data_uvm_cox)
GSE81089_uvm$variable <- c('GSE81089', 'risk_score', 'age', 'stageII', 'stageIII', 'stageIV', 'gender')

data_uvm_cox <- rbind(TCGA_uvm, GSE30219_uvm, GSE31210_uvm, GSE81089_uvm, GSE3141_uvm)
colnames(data_uvm_cox) <- c("variable", "coef", "HR", "lower_95.CI", "upper_95.CI", "pvalue")
data_uvm_cox$mean <-  apply(data_uvm_cox[, c('lower_95.CI', 'upper_95.CI')], 1, function(x){mean(x)})
text_score <-data_uvm_cox[, c('variable', 'HR', 'lower_95.CI', 'upper_95.CI', 'pvalue')]
rmeta_score <- data_uvm_cox[, c("lower_95.CI", "mean", "upper_95.CI")]
text_score <- rbind(colnames(text_score), text_score)
rmeta_score <- rbind(NA, rmeta_score)

pdf(file = 'mut_source_uvm.pdf', 8, 8, onefile = FALSE)
forestplot(text_score, rmeta_score, new_page = TRUE,
           hrzl_lines = list("1" = gpar(lty = 1, lwd = 1.5), 
                             "37" = gpar(lty=1, lwd = 1.5), 
                             '2' = gpar(lty = 2, lwd = 1,  columns=1:5), 
                             '15' = gpar(lty = 2, lwd = 1, columns=1:5), 
                             '22' = gpar(lty = 2, lwd = 1, columns=1:5), 
                             '28' = gpar(lty = 2, lwd = 1, columns=1:5), 
                             '35' = gpar(lty = 2, lwd = 1, columns=1:5)),
           
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 1.2),
                            xlab  = gpar(fontfamily = "", cex = 1.2)),
           
           is.summary=c(TRUE, TRUE, rep(FALSE, 12), TRUE, rep(FALSE, 6), TRUE, rep(FALSE, 4),FALSE, TRUE,  rep(FALSE, 6),  TRUE, FALSE),
           clip = c(0.1, 10), boxsize = 0.2, cex = 3,
           xlog = TRUE, xlab="HR (95% CI)", 
           col=fpColors(box="green4",line="green4", summary="green4", hrz_lines = "#444444"))
dev.off()
rm(text_score, rmeta_score, data_uvm_cox)

###mutivariate cox analysis
load(paste(data_dir, 'TCGA', '/cox/', 'TCGA', '_mvm_cox.Rdata', sep = ''))
TCGA_mvm <- rbind(NA, data_mvm_cox)
TCGA_mvm$variable <- c('TCGA', 'risk_score', 'age', 'T2', 'T3', 'T4', 'N1-3', 'M1', 'stageII', 'stageIII', 'stageIV', 'recurrence', 'gender')

load(paste(data_dir, 'GSE30219', '/cox/', 'GSE30219', '_mvm_cox.Rdata', sep = ''))
GSE30219_mvm <- rbind(NA, data_mvm_cox)
GSE30219_mvm$variable <- c('GSE30219', 'risk_score', 'age' , 'T2', 'N1', 'recurrence', 'gender')

load(paste(data_dir, 'GSE31210', '/cox/', 'GSE31210', '_mvm_cox.Rdata', sep = ''))
GSE31210_mvm <- rbind(NA, data_mvm_cox)
GSE31210_mvm$variable <- c('GSE31210', 'risk_score',  'age', 'stageII', 'recurrence', 'gender')

load(paste(data_dir, 'GSE81089', '/cox/', 'GSE81089', '_mvm_cox.Rdata', sep = ''))
GSE81089_mvm <- rbind(NA, data_mvm_cox)
GSE81089_mvm$variable <- c('GSE81089', 'risk_score', 'age', 'stageII', 'stageIII', 'stageIV', 'gender')

data_mvm_cox <- rbind(TCGA_mvm, GSE30219_mvm, GSE31210_mvm, GSE81089_mvm)
colnames(data_mvm_cox) <- c("variable", "coef", "HR", "lower_95.CI", "upper_95.CI", "pvalue")
data_mvm_cox$mean <-  apply(data_mvm_cox[, c('lower_95.CI', 'upper_95.CI')], 1, function(x){mean(x)})
text_score <-data_mvm_cox[, c('variable', 'HR', 'lower_95.CI', 'upper_95.CI', 'pvalue')]
rmeta_score <- data_mvm_cox[, c("lower_95.CI", "mean", "upper_95.CI")]
text_score <- rbind(colnames(text_score), text_score)
rmeta_score <- rbind(NA, rmeta_score)

pdf(file = 'mut_source_mvm.pdf', 8, 8, onefile = FALSE)
forestplot(text_score, rmeta_score, new_page = TRUE,
           hrzl_lines = list("1" = gpar(lty = 1, lwd = 1.5), 
                             "35" = gpar(lty=1, lwd = 1.5), 
                             '2' = gpar(lty = 2, lwd = 1,  columns=1:5), 
                             '15' = gpar(lty = 2, lwd = 1, columns=1:5), 
                             '22' = gpar(lty = 2, lwd = 1, columns=1:5), 
                             '28' = gpar(lty = 2, lwd = 1, columns=1:5)),
           
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex = 1.2),
                            xlab  = gpar(fontfamily = "", cex = 1.2)),
           
           is.summary=c(TRUE, TRUE, rep(FALSE, 12), TRUE, rep(FALSE, 6), TRUE, rep(FALSE, 4), FALSE, TRUE,  rep(FALSE, 6)),
           clip = c(0.1, 10), boxsize = 0.2, 
           xlog = TRUE, xlab="HR (95% CI)", 
           col = fpColors(box="green4",line="green4", summary="green4", hrz_lines = "#444444"))
dev.off()

