###

library(ggplot2)
library(plyr)
library(gplots)
library(plotly)
library(gmodels)  # the packages that PCA and boxplot need
library(gtools)  # heatmap
library(limma)  #normalizeQuantiles

gnam <- c("ENSG00000019169.9", "ENSG00000122852.10", "ENSG00000168484.8", "ENSG00000185303.11", "ENSG00000203878.7", "ENSG00000231322.1",
          "ENSG00000260695.1", "ENSG00000248608.2", "ENSG00000168878.12", "ENSG00000203878.7", "ENSG00000112175.6", "ENSG00000166961.10")
nam <- "RP11-513N24.1"
nnam <- c("Muscle", "Blood Vessel", "Heart", "Adipose Tissue", "Uterus", "Vagina",  "Breast", "Skin", "Salivary Gland",
          "Brain",  "Adrenal Gland", "Thyroid", "Lung", "Pancreas", "Esophagus", "Stomach", "Colon", "Small Intestine",
          "Prostate",  "Testis", "Nerve", "Blood", "Spleen", "Pituitary", "Ovary", "Liver", "Kidney", "Fallopian Tube",
           "Bladder", "Cervix Uteri")

rt <- read.table(file = "/Users/stead/Desktop/GDC/GTEx/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", header = T, row.names = 1, sep = "\t")
sam_ano <- read.table(file = "/Users/stead/Desktop/GDC/GTEx/GTEx_Data_V6_Annotations_SampleAttributesDS.txt", header = T, row.names = NULL, sep = "\t", quote = "")
rt_ano <-  sam_ano[, c("SAMPID", "SMTS")]
rt_ano[, 1] <- gsub("-", ".", rt_ano[, 1])

rt_g <- data.frame(t(rt[gnam, ])[-1, ])
rt_g_m <- cbind(row.names(rt_g), rt_g)

rt_unit <- merge(rt_ano, rt_g_m, by = 1)
rt_unit <- rt_unit[-c(which(rt_unit[, 2] == "")), ]
#write.table(rt_unit, file = "GATK_12_FPKM.tsv", row.names = T, col.names = T, sep = "\t")
rt_MARCO <- cbind(nam = "MARCO", as.character(rt_unit$SMTS), rt_unit$ENSG00000019169.9)
rt_SFTPA1 <- cbind(nam = "SFTPA1", as.character(rt_unit$SMTS), rt_unit$ENSG00000122852.10)
rt_SFTPC <- cbind(nam = "SFTPC", as.character(rt_unit$SMTS), rt_unit$ENSG00000168484.8)
rt_SFTPA2 <- cbind(nam = "SFTPA2", as.character(rt_unit$SMTS), rt_unit$ENSG00000185303.11)
rt_CHIAP2 <- cbind(nam = "CHIAP2",  as.character(rt_unit$SMTS), rt_unit$ENSG00000203878.7)
rt_RPL13AP17 <- cbind(nam = "RPL13AP17", as.character(rt_unit$SMTS), rt_unit$ENSG00000231322.1)

rt_M <- data.frame(rbind(rt_MARCO, rt_SFTPA1, rt_SFTPC, rt_SFTPA2, rt_CHIAP2, rt_RPL13AP17), stringsAsFactors = F)

#rt_unit <- rt_unit[-c(which(rt_unit[, 3] == 0)), ]

colnames(rt_M) <- c("nam", "group", "exp")

## dimnames <- list(rownames(data), colnames(data))
group <- rt_M[, 2]
rt_M$exp <- log2(as.numeric(rt_M$exp) + 0.001)

### we must define the colors when you have many group
color_group <- c("mediumpurple4", "green", "yellow", "blue", "cyan3", "mediumpurple4", "green", 
                 "yellow", "blue", "cyan3", "mediumpurple4", "green",  "yellow", 
                 "blue", "red", "blue", "yellow", "green", "mediumpurple4", "cyan3", 
                 "blue", "yellow", "green", "mediumpurple4", "cyan3", "blue", "yellow", "yellow",
                "green",  "mediumpurple4")

 
filename <- "GTEx_Boxplot"
tiff(file = paste(filename, ".tiff", sep = ""), width = 2600, height = 3600)
  p = ggplot(data = rt_M, aes(group, rt_M$exp)) +
    facet_wrap(~nam, ncol = 1, nrow = 6, drop = T) +
    xlab("Tissue") +  ylab("gene expression (log2FPKM)")+
    geom_boxplot(color="gray30", fill = group) +
    theme_bw()+
    theme(plot.title = element_text(size=60, hjust = 0.5), 
          strip.text = element_text(size = 50), 
          axis.title.x = element_text(size = 60, color = "black", vjust = 0.5, hjust = 0.5, angle = 1),
          axis.text.x = element_text(size = 60, color = "black", vjust = 0.5, hjust = 0.5, angle = 90),
          axis.title.y = element_text(size = 60, color = "black", vjust = 0.5, hjust = 0.5, angle = 450),
          axis.text.y = element_text(size = 60, color = "black", vjust = 0.5, hjust = 0.5, angle = 1),
          panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
          panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5))
  print(p)
  dev.off()   
