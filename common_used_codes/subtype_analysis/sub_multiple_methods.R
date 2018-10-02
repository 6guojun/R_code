###
###usage: PerSubMain(rt_miRNA_d, rt_sur_m, method_sub = 'hc', rank_sub = 4, data_type = 'miRNA', nrun = 4, seed = 1262118388.71279,  distance = 'pearson')
library(scatterplot3d)
library(plot3D)

source('/Users/stead/Documents/SourceTree/R/common_used_codes/survival_analysis/survival_miner.R')
source('/Users/stead/Documents/SourceTree/R/common_used_codes/survival_analysis/ggsurv_theme/ggsurv_theme_A.R')
source('/Users/stead/Documents/SourceTree/R/common_used_codes/cluster_figure/scatt_3d_PCA.R')

#do heatmap
SubHeat <- function(group_sub, rt_sub_ord, data_type, method_sub, num_k, color_type){
  
  print('group_sub and the expression mat must have the same samples names')
  ha1 = HeatmapAnnotation(group_sub = group_sub, 
                          col = list(group_sub= structure(names = unique(group_sub), color_type[1:length(unique(group_sub))])))
  pdf(file = paste(data_type, method_sub, num_k, 'ht.pdf', sep = '_'), 10, 8)
  ht_list <- Heatmap(rt_sub_ord, name = "center scaled expression",
                     col = colorRamp2(c(min(rt_sub_ord)/3 , (min(rt_sub_ord) + max(rt_sub_ord))/3, max(rt_sub_ord)/2), c(" blue", "white", "red")),
                     top_annotation = ha1, top_annotation_height = unit(0.5, "cm"), 
                     cluster_columns = FALSE, column_dend_reorder = FALSE, 
                     cluster_rows = TRUE, row_dend_reorder = FALSE, 
                     show_row_dend = TRUE, show_column_dend = FALSE,
                     show_row_names = FALSE, show_column_names = FALSE, row_names_side = "right") 
  draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
  dev.off()
  
  print('heatmap finished')
}

###survival analysis
SubSurv <- function(group_sub, rt_sub_ord, rt_sur_m, data_type, method_sub, num_k, color_type){
  int_id <- intersect(colnames(rt_sub_ord), rt_sur_m$sample)
  rt_sur_int <- rt_sur_m[, c('X_OS_IND', 'X_OS')][match(int_id, rt_sur_m$sample, nomatch = 0), ]
  rt_sub_int <- rt_sub_ord[, match(int_id, colnames(rt_sub_ord), nomatch = 0)]
  group_int <- group_sub[match(int_id, names(group_sub), nomatch = 0)]
  
  rt_sub_sur <- data.frame(cbind(colnames(rt_sub_int), rt_sur_int, group_int))
  subtye_name <- paste(data_type, method_sub, num_k, sep = '_')
  colnames(rt_sub_sur) <- c('samples_id', 'OS_Status', 'OS_Time', subtye_name)
  
  DrawSurminer(subtye_name, rt_sub_sur, DatType = "LogType", theme_surv_a, s_cutoff = NULL,  
               palette = color_type[1:length(unique(group_sub))], risk.table = TRUE, fun = NULL)
  print('survival finished')
}

SubPCA <- function(group_sub, rt_sub_ord, data_type, method_sub, num_k, color_type){
  rt_pca <- data.frame(t(rt_sub_ord), stringsAsFactors = FALSE)
  target_genes <- colnames(rt_pca)
  nam <-  paste(data_type, method_sub, num_k, sep = '_')
  sub_color <- as.character(factor(group_sub, labels = color_type[1:length(unique(group_sub))]))
  
  rt_sam <- data.frame(cbind(names(group_sub), group_sub, sub_color))
  colnames(rt_sam) <-  c('samples', 'group', 'color')
  MakPCA(rt_pca, rt_sam = rt_sam, target_genes, nam = nam, phi = 40, pch = 18, len_a = len_a, len_b = len_b,
         len_position = 'bottom', width = NULL, height = NULL)
}


PerSubMain <- function(num_k, rt_sub, rt_sur_m, method_sub = method_sub, rank_sub = rank_sub, data_type = data_type, nrun = 10, seed = 1262118388.71279, 
                       distance = 'pearson'){
  
  mat_sub <- as.matrix(rt_sub)
  
  if(method_sub %in% c('hc', 'km', 'kmdist', 'pam')){
    res_sub = ConsensusClusterPlus(mat_sub, maxK = rank_sub, reps = 100, pItem = 0.8, pFeature = 1, title = paste(data_type, method_sub, sep = '_'),
                                   clusterAlg = method_sub, distance = distance, seed = seed, plot = "pdf")
    print('subtype finished')
    
    sub_type <- as.data.frame(res_sub[[num_k]]['consensusClass'])
    rt_sub_ord <- rt_sub[, match(row.names(sub_type)[order(sub_type$consensusClass)], colnames(rt_sub), nomatch = 0)]
    sub_type_ord <- sub_type[order(sub_type$consensusClass), ]
    names(sub_type_ord) <- row.names(sub_type)[order(sub_type$consensusClass)]
    write.table(sub_type_ord, file = paste(data_type, method_sub, num_k, 'tab.txt', sep = '_'), row.names = TRUE, col.names = TRUE, sep = '\t')
    
    group_sub <- sub_type_ord  
    #do heatmap
    SubHeat(group_sub, rt_sub_ord, data_type, method_sub, num_k, color_type)
    
    #do survival
    SubSurv(group_sub, rt_sub_ord, rt_sur_m, data_type, method_sub, num_k, color_type)
    
    #do pca
    SubPCA(group_sub, rt_sub_ord, data_type, method_sub, num_k, color_type)
    return(group_sub)
    
  } else if (method_sub == 'NMF'){
    sub_er <- nmf(mat_sub, 2:rank_sub, nrun = nrun, seed = seed)
    #do a rank survey
    
    pdf(file = paste(data_type, method_sub, 'nmf_rank_survey.pdf', sep = '_'))
    sub_pl <- plot(sub_er)
    print(sub_pl)
    dev.off()
    
    pdf(file = paste(data_type, method_sub, 'nmf_rank_ht.pdf', sep = '_'))
    cens_pl <- consensusmap(sub_er)
    print(cens_pl)
    dev.off()
    
    for(i in 2: num_k){
      
    res_sub_s <- nmf(mat_sub, i, nrun = nrun, .options = 't')
    sub_type <- predict(res_sub_s)
    rt_sub_ord <- rt_sub[, match(names(sub_type)[order(sub_type)], colnames(rt_sub), nomatch = 0)]
    sub_type_ord <- sub_type[order(sub_type)]
    write.table(sub_type_ord, file = paste(data_type, method_sub, i, 'tab.txt', sep = '_'), row.names = TRUE, col.names = TRUE, sep = '\t')

    group_sub <- sub_type_ord
    
    #do heatmap
    SubHeat(group_sub, rt_sub_ord, data_type, method_sub, i, color_type)
    
    #do survival
    SubSurv(group_sub, rt_sub_ord, rt_sur_m, data_type, method_sub, i, color_type)
    
    #do pca
    SubPCA(group_sub, rt_sub_ord, data_type, method_sub, i, color_type) 
    
    return(group_sub)
    }
  }
  }
