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
  x1 = xc[, 5:10]
  #x2 = xc[, 1:2]
  #apply(x2, 2, table)
  ### can rescale you parameter, when you data have NA or NULL, you can use it
  #xc = covariates(x1, x2, rescale = T, fill.missing = T)
  fc = coxph(y ~ ., data = x1)
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
x <- c(1:6)
coxout(x)


rt_score <- x.c[, c("sampleID", "time", "status", "RFS_time", "RFS_stat", "M", "N", "T","stage", goffm)]
rt_score$`RP11-513N24.11` <- rt_score$`RP11-513N24.11`*(-0.013)
rt_score$`RP11-206P5.2` <- rt_score$`RP11-206P5.2`*(-0.066)
rt_score$SFTPB <- rt_score$SFTPB*(-0.029)
rt_score$CHIAP2 <-rt_score$CHIAP2*(-0.028)
rt_score$BMP5 <- rt_score$BMP5*(-0.087)
rt_score$MS4A15 <- rt_score$MS4A15*(-0.029)



score <- c( rt_score$`RP11-206P5.2` + rt_score$SFTPB + rt_score$BMP5 +
              rt_score$`RP11-513N24.11`+ rt_score$CHIAP2+ rt_score$MS4A15)
rt_score_A <- cbind(rt_score, score)
write.table(rt_score_A, file = "risk_score_data.tsv", row.names = T, col.names = T, sep = "\t")


rt_score_A_S <- rt_score_A[order(rt_score_A$score), ]
#G_score <- c(rep("high_risk", length(which(rt_score_A_S$score >= median(rt_score_A_S$score)))),
#             rep("low_risk", length(which(rt_score_A_S$score < median(rt_score_A_S$score)))))  
G_score <- c(rep("low_risk", length(which(rt_score_A_S$score < quantile(rt_score_A_S$score)[4]))), 
             rep("high_risk", length(which(rt_score_A_S$score >= quantile(rt_score_A_S$score)[4]))))  

G_num <- c(which(rt_score_A_S$score < quantile(rt_score_A_S$score)[4]), which(rt_score_A_S$score >= quantile(rt_score_A_S$score)[4]))
rt_score_N <- cbind(rt_score_A_S, G_score, G_num)

rt_score_S1 <- rt_score_N[which(rt_score_N$stage == 1), ]
rt_score_S2 <- rt_score_N[which(rt_score_N$stage == 2), ]
rt_score_N0 <- rt_score_N[which(rt_score_N$N == 1), ]
rt_score_N1 <- rt_score_N[which(rt_score_N$N == 2), ]
rt_score_M0 <- rt_score_N[which(rt_score_N$M == 1), ]
rt_score_M1 <- rt_score_N[which(rt_score_N$M == 2), ]
filenames <-"M1"
count_sur <- function(x){
  a <- as.character(x[, "G_score"])
  diff <- survdiff(Surv(time, status) ~a, data = x)
  pValue <- 1-pchisq(diff$chisq,df=1)
  pValue <- round(pValue, 5)
  fit <- survfit(Surv(time, status) ~ a, data = x)
  summary(fit)   
  pdf(file = paste(filenames,  ".pdf", sep = ""))
  par(mar = c(5, 5, 5, 3))
  plot(fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
       main = filenames)
  legend(0, 0.17, c("highrisk", "lowrisk"), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
  legend("topright", paste("p-value:", pValue, sep=" "), cex = 1.5)
  dev.off() 
}
count_sur(rt_score_M1)

rt_score_N$status <- gsub(1, "D", rt_score_N$status)
rt_score_N$status <- gsub("0", "S", rt_score_N$status)

color_G <- factor(rt_score_N$G_score, labels = c("red", "blue"))
color_sat <- factor(rt_score_N$status, labels = c("red", "blue"))


tiff(filename = "point_risk_line.tiff", width = 1800, height = 600)
ggplot(rt_score_N, aes(G_num, score)) + geom_point(colour = color_G, size = 5) +
  #geom_vline(xintercept = quantile(G_num)[4], type =2) +
  ylab("risk score")+
  theme_bw()+ 
  theme(plot.title = element_text(size=60, face="bold", hjust = 0.5), 
        strip.text = element_text(size = 60),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 50, color = "black",  vjust = 0.5, hjust = 0.5, angle = 1),
        panel.grid.major = element_blank(), 
        #panel.border = element_rect(color='black', size = 1.5),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=1, colour = "black"))#, axis.ticks = element_line(color = "black", size = 1.5),
dev.off() 



### ggplot make heatmap
#data_heat <- rt_score_N[,  c("sampleID")]
#data_heat_M <- merge(data_heat, rt_score, by = 1)
#data_heat_M_M <- data_heat_M[, c("x", goffm)]
#data_heat_GE_Z <- apply(data_heat_GE, 2, FUN = function(x){(x-median(x))})#/sd(x)
#data.m<- melt(data_heat_M_M)
#p <- ggplot(data.m, aes(x, variable)) + geom_tile(aes(fill = value), colour = "green") + 
#  scale_fill_gradient(low = "blue", high = "red")

### heatmap2 make heatmap
data_heat <- rt_score_N[,  c("sampleID", "G_score", "G_num")]
data_heat_M <- merge(data_heat, x.c, by = 1)
data_heat_L <- data_heat_M[, c("G_score", "G_num", goffm)]
data_heat_L$G_score <- as.character(data_heat_L$G_score)
data_heat_L[, 1][c(data_heat_L[, 1] == "high_risk")] <- 'red'
data_heat_L[, 1][c(data_heat_L[, 1] == "low_risk")] <- 'blue'
data_heat_L <- data_heat_L[order(data_heat_L$G_num), ]
heatmap.cols <- as.character(data_heat_L[, 1])

data_heat_L[, -c(1:2)] <- apply(data_heat_L[, -c(1, 2)], 2, FUN = function(x){(x-median(x))/sd(x)})

data_heat_t <- t(data_heat_L[, -c(1:2)])
#data_heat_t_h <- t(data_heat_L[, -1][which(data_heat_L[, 1] == "red"), ])
#data_heat_t_l <- t(data_heat_L[, -1][which(data_heat_L[, 1] == "blue"), ])

my_palette <- colorRampPalette(c("green", "black", "red"))(n = 1000)
tiff(file = "risk_score_heatmap.tiff",width = 1800, height = 900)
par(mar = c(23, 3, 2, 1))
heatmap.2(data_heat_t, col = colorpanel(99, "blue", "black", "red"), dendrogram = "column",
          key=T, symkey=F, key.xlab = NA, key.title = "Normalied expression", keysize = 1, Colv = F, Rowv = T, trace = "none", density.info = "none",
          labCol = NA, ColSideColors = heatmap.cols, cexRow = 2, colCol = heatmap.cols, scale = "none",  margin = c(14,10), srtCol = 90) 
legend("left", legend = c('high_risk', 'low_risk'), pch = 15, col = c('red', 'blue'), cex = 4, text.col = "black") 
dev.off()

