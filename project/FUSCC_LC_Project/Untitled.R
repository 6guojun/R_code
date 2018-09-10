setwd("/Users/stead/Desktop/20180105_volcano")
rt <- read.table(file = "EXP_RM.txt", header = TRUE, sep = "\t")
rt_c <- rt[, -1]
row.names(rt_c) <- row.names(rt)
rt_d_l <- data.frame(apply(rt_c, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
rt_d_l$log2FC <- rt_d_l[, 2] - rt_d_l[, 1]


DrawSig <- function(num, all_Pval_FC){
  Data_Gnam <- all_Pval_FC[num, ]
  if(Data_Gnam[, "log2FC"] > 1.5)
    return("Up") 
  if(Data_Gnam[, "log2FC"] < -1.5)
    return("Down")
  if(abs(Data_Gnam[, "log2FC"]) <= 1.5)
    return("NoDiff")
}   

GetVolData <- function(all_Pval_FC){
  Sig_L <- lapply(c(1:length(all_Pval_FC[, 1])), DrawSig, all_Pval_FC)
  Sig_C <- c(do.call(rbind, Sig_L))
  all_Pval_FC$Sig <- Sig_C
  return(all_Pval_FC)
} 

all_Pval_FC <- GetVolData(rt_d_l)
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
pdf(file=  "Volcano.pdf")
volcano <- ggplot(all_Pval_FC, aes(x = sample1_tpm, y = sample2_tpm)) +
  geom_point(aes(color = Sig), size = I(2)) +
  scale_color_manual(values =c("blue", "grey", "red")) +
  labs(x = "sample1", y = "sample2") +
  theme_D
print(volcano)
dev.off()



