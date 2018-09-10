setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
gene_list <- c("EDNRB", "S1PR1", "ADRB2", "GRIA1", "CHRNA2", "CHRM1", "ADRA1A", "RXFP1", "CALCRL", "CHRM2", "GRIK4", "ADRB1", "NMUR1", "GHR", "RXFP2", "PTH1R", "C5AR1", "LEPR", "VIPR1", "P2RY1", "PTGER4", "AGTR1", "FPR2", "AVPR2", "ADRA1D", "SSTR1", "S1PR4", "AGTR2", "PTGFR", "ADCYAP1R1", "GRID1", "P2RX2", "P2RY14", "THRA", "P2RX6", "PTGIR", "GRM3", "FPR1", "THRB", "S1PR5", "CNR1", "P2RX7", "NPY5R", "CHRNA4", "P2RX1", "CHRM4", "CHRM3", "PTAFR", "ADRA2C", "NPY1R", "P2RY13", "CTSG", "TBXA2R", "GLP2R", "F2RL3", "CYSLTR1", "CYSLTR2", "NR3C1", "EDNRA", "GRID2", "C3AR1", "APLNR", "TACR1", "NMBR", "SCTR", "NMUR2", "DRD1", "ADRA1B", "TACR3", "GLP1R", "PTGER3", "GABRB2", "PRL", "TSPO", "PTGDR", "CHRNA10", "NPY4R", "PLG", "PTGER2", "GHRHR", "TACR2", "HTR1F", "GRIA3", "LPAR3", "GRIN3B", "FPR3", "ADORA2A", "HCRTR2", "F2R", "S1PR2", "NPFFR1", "GRM1", "PTGER1", "HCRTR1", "GABRB3", "LPAR6", "HRH2", "RXFP4", "GRIN2B", "MAS1", "GRPR", "GABRA2", "GZMA", "HTR7", "GRIN3A", "HTR2A", "HRH4", "ADRB3", "GRM7", "HRH1", "GRIN2A", "LPAR1", "GRIK1", "P2RY2", "GH1", "VIPR2", "CHRNA3", "ADRA2B", "HTR1B", "GRM8", "DRD2", "DRD5", "P2RY8", "AVPR1A", "GABRP", "CHRNA7", "OXTR", "CRHR2", "GRM6", "GRIK3", "OPRL1", "CCKBR", "GRIN2C", "CHRNE", "ADORA1", "ADRA2A", "GRM5", "P2RY11", "PTH2R", "HTR4", "GPR156", "CHRM5", "GABBR1", "OPRK1", "PRLR", "S1PR3", "GABRE", "GRIA2", "P2RY10", "CALCR", "LPAR4", "GPR83", "GRIK5", "GABBR2", "BDKRB2", "LTB4R", "NPFFR2", "DRD4", "CHRNA1", "CHRNB2", "GIPR", "CNR2", "MLNR", "BDKRB1", "GLRB", "PARD3", "F2RL2", "UTS2R", "HTR2B", "NTSR1", "GRIA4", "HTR6", "GABRG3", "GNRHR", "P2RX3", "GABRR1", "P2RX4", "SSTR3", "P2RY4", "LTB4R2", "GRM2", "LEP", "LHCGR", "SSTR2", "CGA", "CHRNA9", "CRHR1", "CHRND", "GABRR2", "OPRD1", "TRPV1", "ADORA2B", "MC5R", "NPBWR1", "CHRNG", "GABRD", "GALR2", "GABRQ", "GRM4", "AVPR1B", "HTR1D", "PRSS1", "F2RL1", "F2", "GABRA3", "PRSS2", "GRIK2", "P2RX5", "MC4R", "LHB", "GCGR", "GRIN2D", "MC1R", "GLRA3", "MTNR1A", "MCHR1", "CHRNA6", "PRSS3", "GPR35", "TSHR", "CHRNB1", "GRIN1", "LPAR2", "CHRNB4", "CHRNA5", "P2RY6", "KISS1R", 
               "TGFBR2", "PLA2G4F", "FGF10", "RASGRP4", "FGFR4", "CACNA2D2", "BDNF", "RASGRF1", "CACNA1S", "FGF2", "PTPN5", "FGFR2", "ARRB1", "RRAS", "CACNB4", "MAP3K3", "RAP1A", "MRAS", "PPP5D1", "DUSP1", "PPP3CC", "ARRB2", "JUND", "RASGRP2", "MAPK10", "PDGFB", "DUSP3", "NTRK2", "IL1A", "AKT3", "GADD45B", "CRK", "FGF18", "NFATC3", "NFATC1", "FGF7", "NR4A1", "RPS6KA2", "ZAK", "TGFB2", "MAP3K8", "FGF14", "DUSP8", "CACNG4", "PLA2G4C", "FOS", "JUN", "RAPGEF2", "RPS6KA1", "SRF", "CACNG1", "IL1B", "IL1R1", "CACNA2D3", "CACNG6", "PRKCB", "PDGFA", "MKNK2", "MAP2K5", "MECOM", "RPS6KA5", "FLNA", "FGFR3", "MAP3K4", "CACNA1C", "MAPK3", "PRKACA", "TGFB1", "MAP4K2", "MEF2C", "FAS", "FGFR1", "HSPA2", "NTF3", "TRAF6", "MAPK11", "TAOK2", "CACNA1D", "MAPT", "TAOK3", "DUSP7", "PTPRR", "RAP1B", "FLNC", "HSPA1L", "MAP3K14", "DUSP16", "RAC2", "NTF4", "STK4", "PDGFRA", "CD14", "NTRK1", "TNF", "CACNA2D1", "FGF9", "TAB2", "PPP3CA", "MAP3K11", "MAPK1", "CACNA1H", "FGF12", "DUSP6", "MAP2K3", "GRB2", "RRAS2", "MAP3K6", "CDC42", "SOS2", "GNG12", "LAMTOR3", "ECSIT", "CACNB2", "PRKACB", "MAPKAPK3", "PDGFRB", "PPP3CB", "MKNK1", "MAP2K4", "FGF1", "NFKB1", "RPS6KA6", "PLA2G4E", "MAPK7", "MAP3K5", "MYC", "TAB1", "AKT1", "RASGRP3", "FGF22", "MAPK8IP1", "MAP3K12", "PAK2", "NGF", "TGFBR1", "FGF17", "CACNA1F", "EGFR", "TGFB3", "CRKL", "MAP2K7", "ATF4", "MAP4K1", "CACNA2D4", "FASLG", "TNFRSF1A", "DUSP2", "PPM1B", "PPM1A", "JMJD7-PLA2G4B", "MAP2K1", "PRKCA", "RAF1", "RELA", "CACNA1A", "BRAF", "RPS6KA3", "IKBKG", "CACNG8", "RASGRF2", "IL1R2", "FGF16", "RASA2", "FGF8", "RASA1", "PLA2G4B", "MAPK9", "MAX", "GADD45G", "TAOK1", "SOS1", "MAP4K4", "FLNB", "RASGRP1", "NFKB2", "MAP3K1", "FGF13", "RELB", "CACNB1", "CACNA1G", "HRAS", "PPP3R1", "MAPK8", "MAP3K2", "DUSP5", "CHUK", "DUSP10", "MAPKAPK2", "MAPK14", "ELK4", "HSPA8", "NRAS", "GADD45A", "CACNA1B", "MAPK12", "MAP2K2", "GNA12", "CDC25B", "HSPA1A", "NLK", "MAPK8IP2", "HSPA6", "AKT2", "PLA2G4D", "PPP5C", "FGF20", "RPS6KA4", "IKBKB", "TP53", "MAP2K6", "DDIT3", "PTPN7", "MAP4K3", "HSPB1", "CACNA1I", "RAC1", "MAP3K7", "NF1", "HSPA1B", "DAXX", "STK3", "ATF2", "PRKCG", "FGF5", "MAPK8IP3", "MAPKAPK5", "STMN1", "KRAS", "MAP3K13", "ELK1", "CACNA1E", "EGF", "DUSP4", "MAPK13", "PLA2G4A", "CACNB3", "FGF19", "DUSP9", "TRAF2", "CASP3", "FGF11", "RAC3", "PAK1")

eg <- bitr(gene_list, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
  
rt_batch <- read.table(file = "batch_files.txt", header = TRUE, row.names = NULL, sep = "\t", stringsAsFactors = FALSE)
rt <- read.table(file = "FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
rt_d <- rt[-which(apply(rt[, -c(405:413)], 1, median) == 0), ]
rt_l <- data.frame(apply(rt_d, 2, function(x){log2(x + 0.01)}), stringsAsFactors = FALSE)

colnames(rt_l) <- matrix(unlist(strsplit(colnames(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 2]
rt_batch_m <- rt_batch[match(colnames(rt_l), rt_batch$samples, nomatch = NA), ]

colnames(rt_l) <- paste(matrix(unlist(strsplit(colnames(rt_l), "[_]")), ncol = 4, byrow = TRUE)[, 4], rt_batch_m$batch_s, sep = "_")


num_pos_A <- grep("SampleA",colnames(rt_l))
names(num_pos_A) <- rep("SampleA", each = length(num_pos_A))
num_pos_B <- grep("SampleB",colnames(rt_l))
names(num_pos_B) <- rep("SampleB", each = length(num_pos_B))

num_pos_T <- grep("T", colnames(rt_l))
names(num_pos_T) <- rep("Tumor", each = length(num_pos_T))
num_pos_N <- grep("N", colnames(rt_l))
names(num_pos_N) <- rep("Normal", each = length(num_pos_N))

rt_samp <- data.frame(cbind(colnames(rt_l), names(sort(c(num_pos_A, num_pos_B, num_pos_T, num_pos_N), decreasing = FALSE))), 
                      stringsAsFactors = FALSE)
colnames(rt_samp) <- c("samples", "group")

rt_sam_col <- rt_samp
rt_sam_col$color <- NA
rt_sam_col$color[grep("Tumor", rt_samp$group)] <- "red"
rt_sam_col$color[grep("Normal", rt_samp$group)] <- "blue"
rt_sam_col$color[grep("SampleA", rt_samp$group)] <- "purple"
rt_sam_col$color[grep("SampleB", rt_samp$group)] <- "green"
rt_sam_col$batch_s <- rt_batch_m$batch_s
colnames(rt_sam_col) <- c("samples", "group", "color", "batch_s")

row.names(rt_l) <- matrix(unlist(strsplit(row.names(rt_l), "[.]")), ncol = 2, byrow = TRUE)[, 1]
rt_l <- rt_l[eg$ENSEMBL, ]
### make a boxplot of all samples form gene level
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_E.R")
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
MakBoxPlot(rt_l, rt_sam_col, theme_E, type = "tiff", width = 3600, height = 1024)


rt_l_d <- rt_l[, -c(405:413)][-grep("NA", row.names(rt_l)), ]
rt_sam_col_d <- rt_sam_col[-c(405:413), ]

### do hca
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/HCA.R")
MakHCA(rt_l_d, rt_sam_col = rt_sam_col_d, "LC_batch1", type = "tiff", method = 'ward.D', cex = 1, width = 3600, height = 716)

### make a pca plot of all samples 
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R")
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_D.R")
MakPCA(rt_l_d, rt_sam = rt_sam_col_d, ThemePCA = theme_D, type = "tiff", width = 1024, height = 512)

data.a <- t(as.matrix(rt_l_d))  
data.pca <- fast.prcomp(data.a, scale = T)  # do PCA
a <- summary(data.pca) 
tmp <- a$importance  # a include 4 sections which contain importance 
pro1 <- as.numeric(sprintf("%.3f", tmp[2,1]))*100 
pro2 <- as.numeric(sprintf("%.3f", tmp[2,2]))*100 # fetch the proportion of PC1 and PC2
pc <- as.data.frame(a$x)  # convert to data.frame
num_pos <- match(row.names(pc), rt_sam_col_d$samples, nomatch = NA)
pc$group <-  rt_sam_col_d$group[num_pos]
pc$color <-  rt_sam_col_d$color[num_pos]

xlab <- paste("PC1(", pro1, "%)", sep = "")
ylab <- paste("PC2(", pro2, "%)", sep = "")

tiff(file = "FUSCC_LC_PCA_text_batch1.tiff", width = 1024, height = 512)
PCA <- ggplot(pc, aes(PC1, PC2, color = group)) +  scale_colour_manual(values = unique(pc$color)) + 
  geom_point(size = 5) + geom_text(label = paste(row.names(pc)), colour="black", size = 4) +
  labs(x = xlab, y = ylab, title = "PCA") + theme_D
print(PCA)
dev.off()
ggplotly(PCA)



