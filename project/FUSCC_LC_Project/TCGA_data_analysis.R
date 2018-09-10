###TCGA data analysis 
###heatmap, pca
###GSEA related pathway


###sorting data
rt_tcga <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/TCGA_LUAD_FPKM.txt", header = TRUE, row.names = 1, sep = "\t", 
                      stringsAsFactors = FALSE)

sam_G <- matrix(unlist(strsplit(colnames(rt_tcga), "[.]")), ncol = 7, byrow = TRUE)[, 4]

tcga_pos_T1 <- grep("01", sam_G)
names(tcga_pos_T1) <- rep("TCGA_T", each = length(tcga_pos_T1))
tcga_pos_T2 <- grep("02", sam_G)
names(tcga_pos_T2) <- rep("TCGA_T", each = length(tcga_pos_T2))
tcga_pos_N <- grep("11", sam_G)
names(tcga_pos_N) <- rep("TCGA_N", each = length(tcga_pos_N))
rt_samp_T <- data.frame(cbind(colnames(rt_tcga), names(sort(c(tcga_pos_T1, tcga_pos_T2, tcga_pos_N), decreasing = FALSE))), 
                        stringsAsFactors = FALSE)
colnames(rt_samp_T) <- c("samples", "group")
rt_tcga_d <- rt_tcga[-which(apply(rt_tcga, 1, median) == 0), ]
rt_tcga_d_l <- data.frame(apply(rt_tcga_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)
#row.names(rt_tcga_d_l) <- matrix(unlist(strsplit(row.names(rt_tcga_d_l), "[.]")), ncol = 2, byrow = TRUE)[, 1]
###get GSEA pathway related gene list hsa04080
gene_list <- c("ADRB2", "S1PR1", "EDNRB", "GRIA1", "RXFP1", "CALCRL", "GRIK4", "NMUR1", "PTGER4", "VIPR1", "S1PR4", "ADRB1", 
               "CHRM1", "C5AR1", "ADRA1A", "PTH1R", "P2RY1", "GHR", "CHRNA2", "THRA", "ADCYAP1R1", "GRID1", "PTGIR", "RXFP2", "THRB", 
               "P2RY14", "LEPR", "FPR2", "CHRM2", "FPR1", "AVPR2", "P2RX7", "F2RL3", "SSTR1", "CNR1", "P2RX1", "TBXA2R", "APLNR", "PTAFR", 
               "P2RY13", "AGTR1", "ADRA1D", "CTSG", "GRM3", "NR3C1", "P2RX6", "EDNRA", "NMBR", "CHRNA10", "S1PR5", "CYSLTR2", "AGTR2", "TSPO", 
               "C3AR1", "P2RX2", "CHRNA4", "CYSLTR1", "PTGFR", "GRIN3B", "CHRM3", "CHRM4", "PTGDR", "S1PR2", "ADORA2A", "F2R", "MAS1", "TACR1", 
               "ADRA1B", "FPR3", "NMUR2", "ADRA2C", "NPY5R", "DRD1", "LPAR6", "GRID2", "HRH4", "PTGER2", "HCRTR1", "NPY1R", "LEP", "SCTR", "BDKRB1", 
               "PTGER1", "GRIA3", "HRH1", "GRM8", "OXTR", "PTGER3", "GRM1", "GH1", "PLG", "GHRHR", "GLP2R", "TACR2", "NPY4R", "NPFFR1", "HCRTR2", "GRIN3A", 
               "RXFP4", "HRH2", "TACR3", "HTR2A", "HTR1F", "GRIN2B", "LPAR3", "DRD5", "VIPR2", "GLP1R", "PRL", "GRM7", "GRM5", "GZMA", "P2RY8", "GRIA2", "GABRB2", 
               "GRPR", "DRD2", "ADRB3", "GRIA4", "HTR7", "GABRA2", "CHRNA3", "HTR1B", "ADRA2B", "CGA", "GRM4", "P2RY2", "F2", "P2RY4", "GRIK1", "AVPR1A", "CCKBR", "PRSS1", "CNR2", 
               "CHRNE", "GRIK3", "HTR6", "NTSR1", "CHRNB2", "MC5R", "LPAR1", "CHRNA7", "PRSS2", "GPR156", "CHRNA1", "UTS2R", "HTR4", "GABRR1", "GRM6", "CHRND", "CRHR2", "GNRHR", 
               "OPRL1", "GPR83", "GALR2", "CHRM5", "GABRP", "CHRNG", "CALCR", "TSHR", "P2RY11", "PTH2R", "F2RL2", "S1PR3", "GCGR", "GABRB3", "GRIN2A", "GRIN2C", "GRM2", "GABBR1", "CRHR1", "GABRG3", 
               "ADRA2A", "DRD4", "HTR2B", "GRIK5", "P2RY10", "LPAR4", "BDKRB2", "NPFFR2", "MLNR", "AVPR1B", "GABRQ", "GLRA3", "CHRNA9", "PRLR", "GIPR", "GABBR2", "OPRK1", "PRSS3", "ADORA1", "NPBWR1", 
               "LTB4R", "GABRE", "OPRD1", "PARD3", "LHCGR", "GABRA3", "P2RX5", "SSTR2", "P2RX4", "SSTR3", "TRPV1", "P2RX3", "LTB4R2", "MCHR1", "GABRD", "GABRR2", "GLRB", "GRIK2", "MC4R", "KISS1R", 
               "CHRNA5", "HTR1D", "GRIN2D", "LHB", "MTNR1A", "ADORA2B", "GRIN1", "CHRNB4", "F2RL1", "CHRNA6", "GPR35", "MC1R", "CHRNB1", "P2RY6", "LPAR2")
eg <- bitr(gene_list, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
rt_tcga_d_l_04080 <- rt_tcga_d_l[eg$ENSEMBL, ]
rt_tcga_d_l_04080 <- rt_tcga_d_l_04080[-grep("NA", row.names(rt_tcga_d_l_04080)), ]



### do hca
rt_sam_col <- rt_samp_T
rt_sam_col$color <- c(rep("green", each = 59), rep("purple", each = 535))
cols <- rt_sam_col$color

rt_hclust <- hclust(dist(t(rt_tcga_d_l_04080)), method = "ward.D2")
rt_hclust_dend <- as.dendrogram(rt_hclust)
new_order <- order.dendrogram(rt_hclust_dend)
labels(rt_hclust_dend) <- colnames(rt_tcga_d_l_04080)[new_order]
labels_colors(rt_hclust_dend) <- rt_sam_col$color[new_order]


tiff(paste("hsa04080", "_HCA.tiff", sep = ""), width = 4096, height = 2048)
par(mar = c(30, 4, 4, 2), cex = 2)
plot(rt_hclust_dend)
colored_bars(cols, rt_hclust_dend, y_scale = 20, y_shift = -100)
legend("topright", legend = unique(rt_sam_col$group), fill = unique(rt_sam_col$color), 
       text.col = unique(rt_sam_col$color), cex = 3)
dev.off()

###do heatmap
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/heatmap.R")
MakHeatmap(rt_tcga_d_l_04080, rt_sam = rt_samp_T, "TCGA_T", "TCGA_N", type = "tiff", width = 1500, height = 1024)

### PCA analysis
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_tcga_d_l_04080, rt_sam = rt_sam_col, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)

###TCGA DEG analysis
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/DEG_T_test.R")
DEGPvalFC <- DEGTest("TCGA_T", "TCGA_N", rt_tcga_d_l, rt_samp_T, var.equal = FALSE, paired = FALSE, meas = "DEGPvalFC")
DEG_S <- DEGPvalFC[order(DEGPvalFC$log2FC), ]
gene_list <- c(row.names(DEG_S)[1 :50], row.names(DEG_S)[8190: 8240])
rt_tcga_d_l_t <- rt_tcga_d_l[gene_list, ]


###do heatmap
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/heatmap.R")
MakHeatmap(rt_tcga_d_l_t, rt_sam = rt_samp_T, "TCGA_T", "TCGA_N", type = "tiff", width = 1500, height = 1024)


### PCA analysis
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_tcga_d_l_t, rt_sam = rt_sam_col, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)

write.table(gene_list, file = "TCGA_LUAD_Top100_gene.txt", row.names = FALSE)
