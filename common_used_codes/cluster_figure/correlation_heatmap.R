###make a correlation heatmap
###you should prepare a table which contain row.names and colnames
###colnames is the elements which you want to perform correlation analysis
###for example: row.names can be all samples value, colnames can be gene names
###example data: mydata <- mtcars[, c(1,3,4,5,6,7)]
###usage: mak_cor_heatmap(mydata, fnam, width = 6, height = 6, cor_text_size = 8, theme_cor_heatmap)
### fnam is the file name of output, theme_cor_heatmap is the ggplot theme which you want to add
###Email: shangjun@163.com
###20180413

library(reshape2)
library(ggplot2)

mak_cor_heatmap <- function(mydata, fnam, width = 6, height = 6, cor_text_size = 8, theme_cor_heatmap){
 
  #get matrix of logicals the same size of a given matrix with entries TRUE in the lower or upper triangle.
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  #Reorder the correlation matrix
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  cormat <- round(cor(mydata), 2)
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  ###Create a ggheatmap
  pdf(file = paste(fnam, ".pdf", sep = ""), width, height)
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "correlation") +
    geom_text(aes(Var2, Var1, label = value), color = "black",  size = cor_text_size) +
    coord_fixed() + theme_cor_heatmap
  print(ggheatmap)
  dev.off()
}







  
 


  
