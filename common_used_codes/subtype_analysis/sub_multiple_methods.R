###
###usage: PerSubMain(rt_miRNA_d, rt_sur_m, method_sub = 'hc', rank_sub = 4, data_type = 'miRNA', nrun = 4, seed = 1262118388.71279,  distance = 'pearson')


PerSubMain <- function(rt_sub, rt_sur_m, method_sub = method_sub, rank_sub = rank_sub, data_type = data_type, nrun = 10, seed = 1262118388.71279, 
                       distance = 'pearson'){
  #
  color_type <- c("firebrick2", "plum3", "green4", "sandybrown", "gold3", "orange1", "yellow2", "peru")
  mat_sub <- as.matrix(rt_sub)
  if(method_sub %in% c('hc', 'km', 'kmdist', 'pam')){
    res_sub = ConsensusClusterPlus(mat_sub, maxK = rank_sub, reps = 100, pItem = 0.8, pFeature = 1, title = paste(data_type, method_sub, sep = '_'),
                                   clusterAlg = method_sub, distance = distance, seed = seed, plot = "pdf")
    print('subtype finished')
    for(i in 2:rank_sub){
      sub_type <- as.data.frame(res_sub[[i]]['consensusClass'])
      rt_sub_ord <- rt_sub[, match(row.names(sub_type)[order(sub_type$consensusClass)], colnames(rt_sub), nomatch = 0)]
      sub_type_ord <- sub_type[order(sub_type$consensusClass), ]
      names(sub_type_ord) <- row.names(sub_type)[order(sub_type$consensusClass)]
      write.table(sub_type_ord, file = paste(data_type, method_sub, i, 'tab.txt', sep = '_'), row.names = TRUE, col.names = TRUE, sep = '\t')
      
      #do heatmap
      group_sub <- sub_type_ord
      ha1 = HeatmapAnnotation(group_sub = group_sub, 
                              col = list(group_sub= structure(names = unique(group_sub), color_type[1:length(unique(group_sub))])))
      
      ###show value distribution
      pdf(file = paste(data_type, method_sub, i, 'ht.pdf', sep = '_'), 10, 8)
      ht_list <- Heatmap(rt_sub_ord, name = "center scaled expression", col = colorRamp2(c(min(mat_sub)/3 , (min(mat_sub) + max(mat_sub))/3, max(mat_sub)/2), c(" blue", "white", "red")),
                         top_annotation = ha1, top_annotation_height = unit(0.5, "cm"), 
                         cluster_columns = FALSE, column_dend_reorder = FALSE, 
                         cluster_rows = TRUE, row_dend_reorder = FALSE, 
                         show_row_dend = FALSE, show_column_dend = FALSE,
                         show_row_names = FALSE, show_column_names = FALSE, row_names_side = "right") 
      draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
      dev.off()
      print('heatmap finished')
      
      ###survival analysis
      rt_sur_ord <- rt_sur_m[, c('X_OS_IND', 'X_OS')][match(colnames(rt_sub_ord), rt_sur_m$sample, nomatch = 0), ]
      rt_sub_sur <- data.frame(cbind(colnames(rt_sub_ord), rt_sur_ord, group_sub))
      colnames(rt_sub_sur) <- c('samples_id', 'status', 'time', 'subtype')
      pdf(paste(data_type, method_sub, i, 'sur.pdf', sep = '_'), onefile = FALSE)
      rt_sur_tmp <<- rt_sub_sur
      fit <- survfit(Surv(time, status) ~ subtype, data = rt_sur_tmp)
      sur_pt <- ggsurvplot(fit, pval = TRUE)
      print(sur_pt)
      dev.off()
      print('survival finished')
    }
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
    
    for(i in 2:rank_sub){
      res_sub_s <- nmf(mat_sub, i, nrun = nrun, .options = 't')
      sub_type <- predict(res_sub_s)
      rt_sub_ord <- rt_sub[, match(names(sub_type)[order(sub_type)], colnames(rt_sub), nomatch = 0)]
      sub_type_ord <- sub_type[order(sub_type)]
      write.table(sub_type_ord, file = paste(data_type, method_sub, i, 'tab.txt', sep = '_'), row.names = TRUE, col.names = TRUE, sep = '\t')
      
      #do heatmap
      group_sub <- sub_type_ord
      ha1 = HeatmapAnnotation(group_sub = group_sub, 
                              col = list(group_sub = structure(names = unique(group_sub), color_type[1:length(unique(group_sub))])))
      pdf(file = paste(data_type, method_sub, i, 'ht.pdf', sep = '_'), 10, 8)
      ht_list <- Heatmap(rt_sub_ord, name = "center scaled expression", col = colorRamp2(c(min(mat_sub) , (min(mat_sub) + max(mat_sub))/2, max(mat_sub)), c(" blue", "white", "red")),
                         top_annotation = ha1, top_annotation_height = unit(0.5, "cm"), 
                         cluster_columns = FALSE, column_dend_reorder = FALSE, 
                         cluster_rows = TRUE, row_dend_reorder = FALSE, 
                         show_row_dend = FALSE, show_column_dend = FALSE,
                         show_row_names = FALSE, show_column_names = FALSE, row_names_side = "right") 
      draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
      dev.off()
      
      #perform survival analysis
      rt_sub_sur = data.frame(cbind(colnames(rt_sub_ord), rt_sur_m[, c('X_OS_IND', 'X_OS')][match(colnames(rt_sub_ord), rt_sur_m$sample, nomatch = 0), ], group_sub))
      colnames(rt_sub_sur) = c('samples_id', 'status', 'time', 'subtype')
      rt_sur_tmp <<- rt_sub_sur
      pdf(paste(data_type, method_sub, i, 'sur.pdf', sep = '_'), onefile = FALSE)
      fit <- survfit(Surv(time, status) ~ subtype, data = rt_sur_tmp)
      sur_pt <- ggsurvplot(fit, pval = TRUE)
      print(sur_pt)
      dev.off()
    }
  }
}

