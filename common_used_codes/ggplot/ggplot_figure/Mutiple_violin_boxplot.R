############make a box or violin plot###########################
###ggplot2 package required
###you should prepare a expression table or expression table and rt_sam
###if your rt_sam is not null you will get violin plot
###if you data only have two group, you will get boxplot with pvalue
###you shoud also prepare a ggplot theme which make the plot perfecter
###UseMeathod: MulVBplot(gnam, rt, rt_sam = NULL, GP1 = GP1, GP2 = GP2, theme_VB, var.equal = FALSE)
###20170907
###JunShang
###E-mail: shangjunv@163.com

library(tidyr)
##############################################################
#setwd("/Users/stead/Documents/SourceTree/R/RNA_seq/DEG_analysis")
#rt <- read.table(file ="WJC_H9_FPKM.txt")
#rt <- rt[apply(rt, 1, mean) > 0, ]
#rt_sam <- read.table(file = "group.txt", header = TRUE, sep = "\t", row.name = NULL, stringsAsFactors = FALSE)
#DEGExp <- DEGTest("H9D0", "H9D2", rt, rt_sam, var.equal = FALSE, paired = FALSE, meas = "DEGExp")
##############################################################



MulVBplot <- function(gnam, rt, rt_sam = NULL, GP1 = GP1, GP2 = GP2, theme_VB, var.equal = FALSE){
  #This functino can be use to make violin plot and boxplot 
  #if your rt_sam is not null you will get violin plot
  #if you data only have two group, you will get boxplot with pvalue

  PreData <- function(gnam, rt, rt_sam = NULL, GP1 = GP1, GP2 = GP2){
    # This function was used to prepare data which was used to draw violin plot
    # your data must be expression table 
    # if your data only have two group, you can choose rt_sam = NULL & !is.null(GP1) & !is.null(GP2)
    # if your data have over two group, you can set GP1 and GP2 as NULL 
    # if your data hvae over two group and you want to choice two group from your data, 
    # you you can choose rt_sam and !is.null(GP1) & !is.null(GP2)
    # the methods:
    #1. rt_vio <- PreData(gnam, rt, rt_sam = rt_sam, GP1 = NULL, GP2 = NULL)
    #2. rt_vio <- PreData(gnam, rt, rt_sam = rt_sam, "H9D0", "H9D2")
    #3. rt_vio <- PreData(gnam, DEGExp, "H9D0", "H9D2")
    
    if(is.null(rt_sam) & !is.null(GP1) & !is.null(GP2)){
      # your data have olny two group
      
      rt_g <- rt[gnam, ]
      rt_vio <- gather(rt_g)
      rt_vio$key[grep(GP1, rt_vio$key)] <- GP1
      rt_vio$key[grep(GP2, rt_vio$key)] <- GP2
      colnames(rt_vio) <- c("group", "expression")
      rt_vio[, "expression"] <- as.numeric(rt_vio[, "expression"])
      print("data with two group")
      return(rt_vio)
    } else if (!is.null(rt_sam) & is.null(GP1) & is.null(GP2)) {
      # your data contain over two group and you want to keep all group in the return
      
      M_position <- match(colnames(rt), rt_sam$samples, nomatch = NA)
      group_nam <- rt_sam$group[M_position]
      rt_c <- rbind(group_nam, rt)
      rt_vio <- data.frame(t(rt_c[c("1", gnam), ]), stringsAsFactors = FALSE)
      colnames(rt_vio) <- c("group", "expression")
      rt_vio[, "expression"] <- as.numeric(rt_vio[, "expression"])
      print("data with over two group")
      return(rt_vio)
    } else if (!is.null(rt_sam)) {
      # your data contain over two group and you want to keep only two group in the return
      
      M_position <- match(colnames(rt), rt_sam$samples, nomatch = NA)
      group_nam <- rt_sam$group[M_position]
      rt_c <- rbind(group_nam, rt)
      rt_vio <- data.frame(t(rt_c[c("1", gnam), ]), stringsAsFactors = FALSE)
      colnames(rt_vio) <- c("group", "expression")
      gp_num <- c(grep(GP1, rt_vio$group), grep(GP2, rt_vio$group))
      rt_vio <- rt_vio[gp_num, ]
      rt_vio[, "expression"] <- as.numeric(rt_vio[, "expression"])
      print("data with over two group and only return two group")
      return(rt_vio)
    } else if (is.null(rt_sam) & is.null(GP1) $ is.null(GP2)) {
      print("rt_sam, GP1 and GP2 can't be NULL at the same time")
    } else {
      print("print check your parameter")
    }
  }
  
  MakTest <- function(rt_test, GP1 = GP1, GP2 = GP2, var.equal = FALSE){
    GPx <- c(grep(GP1, rt_test$group))
    GPy <- c(grep(GP2, rt_test$group))
    x <- as.numeric(rt_test$expression[GPx])
    y <- as.numeric(rt_test$expression[GPy])
    Obj <- t.test(x, y, var.equal = var.equal, paired = FALSE)
    Pval <- round(Obj$p.value, digits = 4)
    return(Pval)
  }
  
  MakVioPlot <- function(rt_vio, gnam, theme_VB){
    pdf(file = paste(gnam, ".pdf"))
    p <- ggplot(rt_vio, aes(x = group, y = expression, fill = group)) + 
      geom_violin(trim=FALSE) + annotate("text", x = 1.5, y = max(rt_vio$expression), label= paste("P = ", Pval, sep = "")) +
      geom_jitter(shape=16, position=position_jitter(0.2)) + labs(title = gnam) +
      theme_VB
    print(p)
    dev.off()
  } 
  
  MakBoxPlot <- function(rt_vio, gnam, theme_VB, GP1 = GP1, GP2 = GP2, var.equal = FALSE){
    Pval <- MakTest(rt_vio, GP1, GP2, var.equal)
    Pval <- round(Pval, 4)
    pdf(file = paste(gnam, ".pdf"))
    p <- ggplot(rt_vio, aes(x = group, y = expression, fill = group)) +
      geom_boxplot() +  annotate("text", x= 1.5, y = max(rt_vio$expression), label= paste("P = ", Pval, sep = "")) +
      geom_jitter(shape=16, position = position_jitter(0.2)) + labs(title = gnam) +
      theme_VB
    print(p)
    dev.off()
  }
  
  
  if(!is.null(rt_sam) & is.null(GP1) & is.null(GP2)){
    rt_vio <- PreData(gnam, rt, rt_sam = rt_sam, GP1 = NULL, GP2 = NULL)
    MakVioPlot(rt_vio, gnam, theme_VB)
  } else if (is.null(rt_sam) & !is.null(GP1) & !is.null(GP2)){
    rt_vio <- PreData(gnam, rt, GP1 = GP1, GP2 = GP2)
    Pval <- MakTest(rt_vio,  GP1 = GP1, GP2 = GP2, var.equal = var.equal)
    MakBoxPlot(rt_vio, gnam, theme_VB, GP1 = GP1, GP2 = GP2, var.equal = var.equal)
  }else if (!is.null(rt_sam) & !is.null(GP1) & !is.null(GP2)) {
    rt_vio <- PreData(gnam, rt, rt_sam = rt_sam, GP1= GP1, GP2 = GP2)
    Pval <- MakTest(rt_vio,  GP1 = GP1, GP2 = GP2, var.equal = var.equal)
    MakVioPlot(rt_vio, gnam, theme_VB) 
  } 
}
