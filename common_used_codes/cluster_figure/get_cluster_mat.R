###get a cluster matrix 
###usage: MakCluterMat(rt_cluster, method = "complete", cluster_dir = c("Vertical", "Horizontal", "all"))
#get cluter order 
library(ggplot2)
library(ggcorrplot)
library(ggdendro)
library(reshape2)

MakCluterMat <- function(rt_cluster, method = "complete", cluster_dir = c("Vertical", "Horizontal", "all")){
  ### dendro_plot_V <- ggdendrogram(data = ord_dendro_V, rotate = FALSE) 
  ### dendro_plot_H <- ggdendrogram(data = ord_dendro_H, rotate = TRUE)
  print("you will get a list contain reorder matrix and ord_dendro_V or y 
        you can get clustr plot by using 
        dendro_plot_V <- ggdendrogram(data = ord_dendro_V, rotate = FALSE) and 
        dendro_plot_H <- ggdendrogram(data = ord_dendro_H, rotate = FALSE)")
  
  ord_dendro_V <- as.dendrogram(hclust(d = dist(x = t(rt_cluster)), method = method))
  ord_dendro_H <- as.dendrogram(hclust(d = dist(x = rt_cluster),  method = method))
  otter_order_V <- order.dendrogram(ord_dendro_V)
  otter_order_H <- order.dendrogram(ord_dendro_H)
  #Re-order heatmap columns to match dendrogram
  if(cluster_dir == "Vertical"){
    print("you will get a list contain reorder matrix and ord_dendro_V which cluster in column, Vertical direction")
    rt_reorder <- rt_cluster[, otter_order_V]
    reorder_list <- list(rt_reorder, ord_dendro_V)
    names(reorder_list) = c("reorder_mat", "ord_dendro_V")
    return(reorder_list)
  } else if(cluster_dir == "Horizontal"){
    print("you will get a list contain reorder matrix and ord_dendro_V which cluter in row, Horizontal direction")
    rt_reorder <- rt_cluster[otter_order_H, ]
    reorder_list <- list(rt_reorder, ord_dendro_H)
    names(reorder_list) = c("reorder_mat", "ord_dendro_H")
    return(reorder_list)
  } else if (cluster_dir == "all") {
    print("you will get a list contain reorder matrix and ord_dendro_V/H which cluster row and column, Vertical and Horizontal direction")
    rt_reorder <- rt_cluster[otter_order_H, otter_order_V]
    reorder_list <- list(rt_reorder, ord_dendro_V, ord_dendro_H)
    names(reorder_list) = c("reorder_mat", "ord_dendro_V", "ord_dendro_H")
    return(reorder_list)
  } else {
    stop("you must set cluster_dir")
  }
}

