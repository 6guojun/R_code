############make a volcano plot################################
###ggplot2 package required
###you should prepare a table contain pvalue and fold change
###three columns should contain in your table: "Padj", "log2FC" and "Sig" 
###fc_v come from log2(fc_v), so fc_v must > 0
###you shoud also prepare a ggplot theme which make the plot perfecter
###UseMeathod: MakVolPlot(AllPvalFC, GP1, GP2, fc_v, theme_A, width, height)
###20170907
###JunShang
###E-mail: shangjunv@163.com


##############test data#########################################
#source("/Users/stead/Documents/SourceTree/R/RNA_seq/DEG_analysis/DEG_Test.R")
#rt <- read.table(file ="WJC_H9_FPKM.txt")
#rt <- rt[apply(rt, 1, mean) > 0, ]
#rt_sam <- read.table(file = "group.txt", header = TRUE, sep = "\t", row.name = NULL, stringsAsFactors = FALSE)
#AllPvalFC <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "AllPvalFC")
#source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_A.R")
################################################################

######################plot a volcano############################
MakVolPlot <- function(all_Pval_FC, GP1, GP2, fc_v, Vtheme, width, height){
  #DrawSig can tag downexpression or overexpression on each gene
  #GetVolData can get a table which can be used to do volcano plot
  #GP1 and GP2 were just used to produce figure name

  DrawSig <- function(Gnam, all_Pval_FC){
    Data_Gnam <- all_Pval_FC[Gnam, ]
    if(Data_Gnam[, "Padj"] < 0.05 & Data_Gnam[, "log2FC"] > fc_v)
      return("Up")
    if(Data_Gnam[, "Padj"] < 0.05 & Data_Gnam[, "log2FC"] < -fc_v)
      return("Down")
    if(Data_Gnam[, "Padj"] > 0.05 | abs(Data_Gnam[, "log2FC"]) <= fc_v)
      return("NoDiff")
  }
  
  GetVolData <- function(all_Pval_FC){
    Sig_L <- lapply(row.names(all_Pval_FC), DrawSig, all_Pval_FC)
    Sig_C <- c(do.call(rbind, Sig_L))
    all_Pval_FC$Sig <- Sig_C
    return(all_Pval_FC)
  }
  
  all_Pval_FC <- GetVolData(all_Pval_FC)

  pdf(file= paste(GP1, "_", GP2, "_", "Volcano", ".pdf", sep = ""), width, height)
  volcano <- ggplot(all_Pval_FC, aes(x = log2FC, y = -1*log10(Padj))) + 
    geom_point(aes(color = Sig), size = I(2)) + 
    scale_color_manual(values =c("blue", "grey", "red")) +
    labs(title= "Volcanoplot", x = "log2(fold change)", y = "-log10(p value)") +
    Vtheme
  print(volcano)
  dev.off()  
}
################################################################
