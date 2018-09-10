library(ComplexHeatmap)
library(limma)
library(circlize)
library(data.table)
library(ggplot2)
library(ggpubr)

setwd('/Users/stead/Desktop/CHIAP2/miRNA/limma_DEGS')

rt_miRNA <- read.table(file = '/Users/stead/Desktop/CHIAP2/miRNA/miRNA_exp_d.txt', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_miRNA_l <- data.frame(apply(rt_miRNA, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
colnames(rt_miRNA_l) <- c('shCtrl_1', "shCtrl_2", "shCtrl_3", 'shCHIAP2_1', 'shCHIAP2_2', 'shCHIAP2_3')  
rt_miRNA_l <- rt_miRNA_l[, c('shCHIAP2_1', 'shCHIAP2_2', 'shCHIAP2_3', 'shCtrl_1', "shCtrl_2", "shCtrl_3")]

rt_sam_mi <- data.frame(cbind(colnames(rt_miRNA_l), c(rep('shCHIAP2', 3), rep('shCtrl', 3))), stringsAsFactors = FALSE)
colnames(rt_sam_mi) <- c('samples_id', 'group')

mi_group = factor(rt_sam_mi$group, levels=c('shCHIAP2', 'shCtrl'))

design_mi = model.matrix(~0 + mi_group)
row.names(design_mi) <- rt_sam_mi$samples_id
colnames(design_mi) <- c('shCHIAP2','shCtrl')
design_mi

###
mi_fit <- lmFit(rt_miRNA_l,  design_mi)
cont.matrix <- makeContrasts(shCHIAP2vsshCtrl = shCHIAP2-shCtrl, levels = design_mi)
mi_fit2 <- contrasts.fit(mi_fit, cont.matrix)
mi_fit3 <- eBayes(mi_fit2)
mir_diff <- topTreat(mi_fit3, number = 1589)
mir_diff2 <- mir_diff[which(abs(mir_diff$logFC) > 0.5 & mir_diff$P.Value < 0.05), ]
write.table(mir_diff2, file = "limma_miRNA_DEGs.txt", col.names = TRUE, row.names = TRUE, sep = "\t")

rt_miRNA_DEG <- rt_miRNA_l[match(row.names(mir_diff2), row.names(rt_miRNA), nomatch = 0), ]

###do a heatmap
mat_scale_mi <- t(apply(rt_miRNA_DEG, 1, scale))
colnames(mat_scale_mi) <- colnames(rt_miRNA_DEG)

pdf("CHIAP2_miRNA_DEG.pdf", 8, 8)
ha = HeatmapAnnotation(type = c(rep("shCHIAP2", 3), rep("shCtrl", 3)), col = list(type = structure(names = c("shCHIAP2", "shCtrl"), c("red", "blue")))) 
                       
HT = Heatmap(mat_scale_mi, name = "expression", col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(4, "mm"), 
        show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 13, fontface = "bold"),
        cluster_columns = TRUE, column_dend_reorder = FALSE, 
        cluster_rows = TRUE, row_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
        show_row_dend = TRUE, show_column_dend = TRUE, show_heatmap_legend = TRUE)
draw(HT, heatmap_legend_side = "right",  annotation_legend_side = "bottom")
dev.off()

###volcano analysis
###make volcano plot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/volcano_plot.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_D.R")
mir_diff_vol <- mir_diff[, c('P.Value', "logFC")]
colnames(mir_diff_vol) <- c('Padj', 'log2FC')
MakVolPlot(mir_diff_vol, "shCHIAP2", "shCtrl", 0.5, theme_D, 8, 8)

###PCA plot
###make a pca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
rt_sam_mi$color <- c(rep("red", 3), rep("blue", 3))
MakPCA(rt_miRNA_DEG, rt_sam = rt_sam_mi, nam = "miRNA", theme_D, type = "pdf", 6.5, 5)

#------------------
rm(mi_fit, mir_diff, mir_diff2, cont.matrix, mir_diff_vol, rt_miRNA, design_mi, 
   rt_miRNA_DEG, rt_miRNA_l, mat_scale_mi, rt_sam_mi, ha, mi_fit, mi_fit2, mi_fit3)
####################################################mRNA analysis####################################################

###mRNA 
setwd('/Users/stead/Desktop/CHIAP2/mRNA/limma_DEGs')
rt_mRNA <- read.table(file = "/Users/stead/Desktop/CHIAP2/mRNA/FPKM_mRNA_sym.txt", header = TRUE, row.names = 1, sep = "\t")
rt_mRNA_d <- rt_mRNA[apply(rt_mRNA, 1, function(x){mean(x) > 0}), ]
rt_mRNA_l <- data.frame(apply(rt_mRNA_d, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
colnames(rt_mRNA_l) <- c('shCtrl_1', "shCtrl_2", "shCtrl_3", 'shCHIAP2_1', 'shCHIAP2_2', 'shCHIAP2_3')
rt_mRNA_l <- rt_mRNA_l[, c('shCHIAP2_1', 'shCHIAP2_2', 'shCHIAP2_3', 'shCtrl_1', "shCtrl_2", "shCtrl_3")]


rt_sam_m <- data.frame(cbind(colnames(rt_mRNA_l), c(rep('shCHIAP2', 3), rep('shCtrl', 3))), stringsAsFactors = FALSE)
colnames(rt_sam_m) <- c('samples_id', 'group')
group_m = factor(rt_sam_m$group, levels=c('shCHIAP2', 'shCtrl'))
design_m = model.matrix(~0 + group_m)
row.names(design_m) <- rt_sam_m$samples_id
colnames(design_m) <- c('shCHIAP2', 'shCtrl')
design_m

###DEGs
###
fit_m <- lmFit(rt_mRNA_l,  design_m)
cont.matrix <- makeContrasts(shCHIAP2vsshCtrl = shCHIAP2-shCtrl, levels = design_m)
fit_m2 <- contrasts.fit(fit_m, cont.matrix)
fit_m3 <- eBayes(fit_m2)
m_diff <- topTreat(fit_m3, number = length(row.names(rt_mRNA_l)))
m_diff2 <- m_diff[which(abs(m_diff$logFC) > 0.5 & m_diff$P.Value < 0.05), ]
write.table(m_diff2, file = "limma_miRNA_DEGs.txt", col.names = TRUE, row.names = TRUE, sep = "\t")

rt_mRNA_DEG <- rt_mRNA_l[match(row.names(m_diff2), row.names(rt_mRNA_l), nomatch = 0), ]


###do a heatmap
mat_scale_m <- t(apply(rt_mRNA_DEG, 1, scale))
colnames(mat_scale_m) <- colnames(rt_mRNA_DEG)

pdf("CHIAP2_mRNA_DEG.pdf", 8, 8)
ha = HeatmapAnnotation(type = c(rep("shCHIAP2", 3), rep("shCtrl", 3)), col = list(type = structure(names = c("shCHIAP2", "shCtrl"), c("blue", "red"))))
Heatmap(mat_scale_m, name = "expression", col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation = ha, top_annotation_height = unit(4, "mm"), 
        show_row_names = FALSE, show_column_names = TRUE,
        cluster_columns = TRUE, column_dend_reorder = FALSE, 
        cluster_rows = TRUE, row_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 25, fontface = "bold"),
        show_row_dend = TRUE, show_column_dend = TRUE) 
dev.off()

###volcano analysis
###make volcano plot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/volcano_plot.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_D.R")
m_diff_vol <- m_diff[, c('P.Value', "logFC")]
colnames(m_diff_vol) <- c('Padj', 'log2FC')
MakVolPlot(m_diff_vol,"shCHIAP2", "shCtrl", 0.5, theme_D, 5, 5)

###PCA plot
###make a pca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
rt_sam_m$color <- c( rep("red", 3), rep("blue", 3))
MakPCA(rt_mRNA_DEG, rt_sam = rt_sam_m, nam = "mRNA", theme_D, type = "pdf", 6.5, 5)

#------------------
rm(fit_m, fit_m2, m_diff, m_diff2, cont.matrix, m_diff_vol, rt_mRNA, design_m, 
   rt_mRNA_DEG, rt_mRNA_l, mat_scale_m, rt_sam_m, ha, fit_m3, rt_mRNA_d)
#######taget predict
####################################################mRNA analysis####################################################
setwd('/Users/stead/Desktop/CHIAP2/target_predict')
rt_GSEA <- read.table(file = "/Users/stead/Desktop/CHIAP2/mRNA/GSEA/old/gsea_kegg/gseapy.gsea.gene_sets.report.csv", header = TRUE, 
                      sep = ",", stringsAsFactors = FALSE)

wnt_genes_list <- rt_GSEA$genes[which(rt_GSEA$Term == 'KEGG_WNT_SIGNALING_PATHWAY')]
wnt_genes <- unlist(strsplit(wnt_genes_list, ","))
write.table(wnt_genes, file = "wnt_genes.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


###do all insect genes heatmap and boxplot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_B.R")


mir_name <- c('let_7f_5p',	'mir_3129_5p',	'mir_3173_3p',	'mir_3614_5p',	'mir_873_3p')
for(mir_n in mir_name){
  setwd(paste('/Users/stead/Desktop/CHIAP2/target_predict/', mir_n, sep = ""))
  mir_tar_Tscan <- read.table(file = "Tscan_target_genes.txt", header = TRUE, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
  mir_tar_miRDb <- read.table(file = "miRDB_target_genes.txt", header = TRUE, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
  mir_wnt_int <- intersect(intersect(mir_tar_Tscan$Target, wnt_genes), mir_tar_miRDb$Detail)
 print(length(mir_wnt_int))

 rt_mir_wnt_int <- rt_mRNA_l[mir_wnt_int, ]
 mat_scale_m <- t(apply(rt_mir_wnt_int, 1, scale))
 colnames(mat_scale_m) <- colnames(rt_mir_wnt_int)
 
 pdf(paste(mir_n, 'heatmap.pdf', sep = '_'), 5, 8)
 ha = HeatmapAnnotation(type = c(rep("shCtrl", 3), rep("shCHIAP2", 3)), col = list(type = structure(names = c("shCtrl", "shCHIAP2"), c("blue", "red"))))
 pt <- Heatmap(mat_scale_m, name = "expression", col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               top_annotation = ha, top_annotation_height = unit(4, "mm"), 
               show_row_names = TRUE, show_column_names = TRUE,
               cluster_columns = TRUE, column_dend_reorder = FALSE, 
               cluster_rows = TRUE, row_dend_reorder = FALSE, column_names_gp = gpar(fontsize = 18, fontface = "bold"), 
               show_row_dend = TRUE, show_column_dend = TRUE) 
 print(pt)
 dev.off()
 
 rt_box_raw <- melt(t(rt_mir_wnt_int))
 group <- rep(c(rep('shCTRl', 3), rep('shCHIAP2', 3)), length(row.names(rt_mir_wnt_int)))
 rt_box <- data.frame(cbind(rt_box_raw, group), stringsAsFactors = FALSE)
 colnames(rt_box) <- c('samples_id', "gene_name", "expression", "group")
 pdf(file = paste(mir_n, 'boxplot.pdf', sep = '_'), 5, 4)
 e_box_raw <- ggplot(rt_box, aes(x = gene_name, y = expression))
 e_box <- e_box_raw + 
   geom_boxplot(aes(fill = group), position = position_dodge(0.9) ) + 
   scale_fill_manual(values = c("red", "blue")) + 
   stat_compare_means(aes(group = group), label = "p.signif") + 
   theme_classic_A + ylab("log2(exp + 1)")
 print(e_box)
 dev.off()
}




