############make a box plot###########################
###ggplot2, tidyr, dplyr package required
###you should prepare expression table and sample table 
###meanwhile samples table sholud contain group which was accord with the name of samples
###colnames of rt_box should contain gnams(one or more) and "group"
###UseMeathod:MakBoxPlot(gnam, "group", rt_box, width = 4, height = 6, theme_VB) 
###20190415
###JunShang
###E-mail: shangjunv@163.com

library(reshape2)
library(ggplot2)

#================================================================
#convert an epxression table to a dataframe which format can be used by ggplot
#
#================================================================
#the format of rt_exp like tcga expression table 
#the rt_sam contain two columns, one of which is sampls_id and the other is group type
#the samples id of rt_sam must have the same order with colnames of rt_exp
#usage: rt_g <- conver_exp_data(rt_exp, rt_sam)
#col_names can used to set colnames of your output data
conver_exp_data <- function(rt_exp, rt_sam, col_names = c("samples_id", "group", "gene_id", "exp")){
  print(col_names)
  rt_exp_t <- data.frame(t(rt_exp), stringsAsFactors = FALSE)
  rt_exp_m <- data.frame(cbind(rt_sam, rt_exp_t))
  rt_g <- melt(rt_exp_m)
  colnames(rt_g) <- col_names
  return(rt_g)
}




#=================================================================
#make a box plot or violin plot with ggplot2
#
#=================================================================
#you should prepare a matrix which at least contain two colums
#the one colums is group and the other is expression value or other value
#the data format based on the boxplot format of ggplot2 
#you will get two table, one of which is boxplot and the other is violin plot
#usage: MakBoxPlot(exp, group, rt_box,width = 4, height = 6, theme_VB)

MakBoxPlot <- function(gnam, gggroup, rt_box_mat, width = 4, height = 6, theme_VB, ggtype = c('boxplot', 'violin'), ggcolor = ggcolor, ggylab = ggylab){
   #group contain all kinds of types such as tumor or normal, female or male, cancer type and so on. 
   #the column of gnam contain all group expression value or other value
   #rt_box must have two columns which contain gnam and group. But, it can also contian samples_id and so no.
  gnam <<- gnam
  gggroup <<- gggroup
  rt_box <- rt_box_mat[, c(gnam, gggroup)]
  rt_box[, gnam] <- as.numeric(rt_box[, gnam])
  colnames(rt_box) <- c(gnam, "group")
  
  if(length(unique(rt_box$group)) > 2){
      print("ANOVA will be used becuase your group > 2")
      res_aov <- aov(rt_box[, gnam] ~ group, data = rt_box)
      Pval <- round(unlist(summary(res_aov))[9], digits = 4)

      } else if (length(unique(rt_box$group)) == 2){
      #do t.test
      print("t.test will be used becuase your group = 2")
      GPx <- rt_box[, gnam][grep(unique(rt_box$group)[1], rt_box$group)]
      GPy <- rt_box[, gnam][grep(unique(rt_box$group)[2], rt_box$group)]
      Obj <- t.test(GPx, GPy, var.equal = FALSE, paired = FALSE)
      Pval <- round(Obj$p.value, digits = 4)
    } else {
      stop("your group < 2")
    }
  
  ###do boxplot
  if(ggtype == "boxplot"){
    pdf(file = paste(gnam, '_', gggroup,"_boxplot_point.pdf", sep = ""), width = width, height = height)
    p <- ggplot(rt_box, aes(x = group, y = rt_box[, gnam], fill = group)) + geom_boxplot() + 
      scale_fill_manual(values= ggcolor) + annotate("text", x = 1.5, y = max(rt_box[, gnam]), label= paste("P = ", Pval, sep = "")) +  
      geom_jitter(shape=16, position = position_jitter(0.2)) +  labs(title = gnam) + xlab(gggroup) + ylab(ggylab) +
      theme_VB
    print(p)
    dev.off()
    
    #no point 
    pdf(file = paste(gnam, '_', gggroup, "_boxplot.pdf", sep = ""), width = width, height = height)
    p <- ggplot(rt_box, aes(x = group, y = rt_box[, gnam], fill = group)) + geom_boxplot() + 
      scale_fill_manual(values= ggcolor) + annotate("text", x = 1.5, y = max(rt_box[, gnam]), label= paste("P = ", Pval, sep = ""))  +
      labs(title = gnam) + xlab(gggroup) + ylab(ggylab) + theme_VB
    print(p)
    dev.off()
    rm(gnam, gggroup, pos = ".GlobalEnv")
  } else if (ggtype == "violin"){
    pdf(file = paste(gnam, '_', gggroup, "_violin.pdf", sep = ""), width = width, height = width)
    p <- ggplot(rt_box, aes(x = group, y = rt_box[, gnam])) +  geom_violin(aes(fill = group), trim = FALSE) + geom_boxplot(width = 0.1) + 
      scale_fill_manual(values = ggcolor)+ annotate("text", x = 1.5, y = max(rt_box[, gnam]), label= paste("P = ", Pval, sep = "")) +
      labs(title = gnam)  + xlab(gggroup) + ylab(ggylab) + theme_VB
    print(p)
    dev.off()
    rm(gnam, gggroup, pos = ".GlobalEnv")
  } else {
    stop('you must set ggtype')
  }
}
