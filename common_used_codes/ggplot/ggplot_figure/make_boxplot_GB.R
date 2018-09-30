source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_no_xtick.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_E.R")

BGbox1 <- function(rt_box, TCGA_GAS_pvaue, gnam){
  gene_list <- paste(unique(rt_box$variable), collapse ="\n")
  pvalue_vec <- c(TCGA_GAS_pvaue$pvalue[c(which(row.names(TCGA_GAS_pvaue) == gnam), grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue)))])
  names(pvalue_vec) <- c(gnam, row.names(TCGA_GAS_pvaue)[grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue))])
  
  pl_sig <- which(pvalue_vec <= 0.05)
  rt_rect <- data.frame(x1 = c(-0.5 + 1*pl_sig), x2 = c(0.5 + 1*pl_sig), y1 = c(rep(-Inf, length(pl_sig))), y2 = c(rep(Inf, length(pl_sig))))
  gnum = c(1.5 + 1:length(pvalue_vec))
  
  BG_plot <- ggplot2.boxplot(data = rt_box, xName = 'variable', yName = 'value', groupName = 'group', outlier.shape = NA, 
                             groupColors=c('green','red'),  xtitle = gnam, coef = 0, 
                             legendPosition = "bottom", ytitle = "log2(FPKM + 1)",  addDot = FALSE, dotSize = 0.3) +
    stat_boxplot(geom = "errorbar", width = 0.75, size = 0.6,  color = "black") + 
    geom_vline(xintercept = 1.5, color = "black", size = 1) + 
    #geom_hline(yintercept = 1, color = "black", size = 1, linetype = "dashed") + 
    #guides(fill = guide_legend(title = gene_list, title.position = "right", ncol = 4, byrow = TRUE)) + 
    geom_vline(xintercept = gnum, color = "black", size = 0.7, linetype = "dashed") + 
    geom_rect(data = rt_rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = 'pink', alpha = 0.2, inherit.aes = FALSE) +
    #stat_compare_means(aes(group = group), method = "t.test", label = "p.signif") + 
    geom_boxplot(aes(fill = group), outlier.shape = NA, coef = 0) + 
    theme_E  
  print(BG_plot)
}

BGbox2 <- function(rt_box, TCGA_GAS_pvaue, gnam){
  gene_list <- paste(unique(rt_box$variable), collapse ="\n")
  pvalue_vec <- c(TCGA_GAS_pvaue$pvalue[c(which(row.names(TCGA_GAS_pvaue) == gnam), grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue)))])
  names(pvalue_vec) <- c(gnam, row.names(TCGA_GAS_pvaue)[grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue))])
  
  pl_sig <- which(pvalue_vec <= 0.05)
  rt_rect <- data.frame(x1 = c(-0.5 + 1*pl_sig), x2 = c(0.5 + 1*pl_sig), y1 = c(rep(-Inf, length(pl_sig))), y2 = c(rep(Inf, length(pl_sig))))
  gnum = c(1.5 + 1:length(pvalue_vec))
  
  BG_plot <- ggplot2.boxplot(data = rt_box, xName = 'variable', yName = 'value', groupName = 'group', outlier.shape = NA, 
                             groupColors=c('green','red'),  xtitle = gnam,
                             legendPosition = "bottom", ytitle = "log2(FPKM + 1)",  addDot = FALSE, dotSize = 0.3) +
    stat_boxplot(geom = "errorbar", width = 0.75, size = 0.6,  color = "black") + 
    #geom_hline(yintercept = 1, color = "black", size = 1, linetype = "dashed") + 
    #guides(fill = guide_legend(title = gene_list, title.position = "right", ncol = 4, byrow = TRUE)) + 
    geom_vline(xintercept = gnum, color = "black", size = 0.5, linetype = "dashed") + 
    #geom_rect(data = rt_rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = 'pink', alpha = 0.2, inherit.aes = FALSE) +
    #stat_compare_means(aes(group = group), method = "t.test", label = "p.signif") + 
    geom_boxplot(aes(fill = group), outlier.shape = NA) + 
    theme_E  
  print(BG_plot)
}

BGbox1 <- function(rt_box, TCGA_GAS_pvaue, gnam){
  gene_list <- paste(unique(rt_box$variable), collapse ="\n")
  pvalue_vec <- c(TCGA_GAS_pvaue$pvalue[c(which(row.names(TCGA_GAS_pvaue) == gnam), grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue)))])
  names(pvalue_vec) <- c(gnam, row.names(TCGA_GAS_pvaue)[grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue))])
  
  pl_sig <- which(pvalue_vec <= 0.05)
  rt_rect <- data.frame(x1 = c(-0.5 + 1*pl_sig), x2 = c(0.5 + 1*pl_sig), y1 = c(rep(-Inf, length(pl_sig))), y2 = c(rep(Inf, length(pl_sig))))
  gnum = c(1.5 + 1:length(pvalue_vec))
  
  BG_plot <-  ggplot(rt_box, aes(x = variable, y = value, fill = group)) 
  BG_plot_box <- BG_plot + 
    stat_boxplot(geom = "errorbar", width = 0.75, size = 0.6,  color = "black") + 
    geom_boxplot(aes(fill = group), outlier.shape = NA, coef = 0) +
    xlab(gnam) + ylab("expression level log2(FPKM + 1)") + ggtitle(gnam) + 
    geom_rect(data = rt_rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = 'pink', alpha = 0.2, inherit.aes = FALSE) +
    geom_hline(yintercept = 1, color = "black", size = 1, linetype = "dashed") + 
    geom_vline(xintercept = 1.5, color = "black", size = 1) + 
    geom_vline(xintercept = gnum, color = "black", size = 0.5, linetype = "dashed") + 
    geom_boxplot(aes(fill = group), outlier.shape = NA, coef = 0) +
    scale_fill_manual(values=c("green", "red")) + 
    theme_E
  print(BG_plot_box)
}

BGbox2 <- function(rt_box, TCGA_GAS_pvaue, gnam){
  gene_list <- paste(unique(rt_box$variable), collapse ="\n")
  pvalue_vec <- c(TCGA_GAS_pvaue$pvalue[c(which(row.names(TCGA_GAS_pvaue) == gnam), grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue)))])
  names(pvalue_vec) <- c(gnam, row.names(TCGA_GAS_pvaue)[grep(paste(gnam, '_', sep = ""), row.names(TCGA_GAS_pvaue))])
  
  pl_sig <- which(pvalue_vec <= 0.05)
  rt_rect <- data.frame(x1 = c(-0.5 + 1*pl_sig), x2 = c(0.5 + 1*pl_sig), y1 = c(rep(-Inf, length(pl_sig))), y2 = c(rep(Inf, length(pl_sig))))
  gnum = c(1.5 + 1:length(pvalue_vec))
  
  BG_plot <-  ggplot(rt_box, aes(x = variable, y = value, fill = group)) 
  BG_plot_box <- BG_plot + 
    stat_boxplot(geom = "errorbar", width = 0.75, size = 0.6,  color = "black") + 
    geom_boxplot(aes(fill = group), outlier.shape = NA, coef = 0) +
    xlab(gnam) + ylab("expression level log2(FPKM + 1)") + ggtitle(gnam) + 
    #geom_rect(data = rt_rect, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = 'pink', alpha = 0.2, inherit.aes = FALSE) +
    #geom_hline(yintercept = 1, color = "black", size = 1, linetype = "dashed") + 
    geom_vline(xintercept = 1.5, color = "black", size = 1) + 
    geom_vline(xintercept = gnum, color = "black", size = 0.5, linetype = "dashed") + 
    geom_boxplot(aes(fill = group), outlier.shape = NA, coef = 0) +
    scale_fill_manual(values=c("green", "red")) + 
    theme_E
  print(BG_plot_box)
}
