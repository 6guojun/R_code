### you can count any gene expression between tumor and normal
### prepare two file RNA-seq data and miRNA data

library(ggplot2)
cnam <- "LUAD"
gnam <- "ENSG00000166961.13"
rt <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnam, "/raw_data/", cnam, "_mRNA_FPKM.txt", sep = ""), header = T, sep = "\t")
splite_data <- function(rt, k){
  ### splite the tumor and normal base on k which is a character 
  
  keepRow <- c(1: nrow(rt))            
  keepCol <- seq(1, ncol(rt)) ## product the data of read per millions
  data <- rt[keepRow, keepCol]
  data <- log2(data + 0.001)
  dimnames <- list(rownames(data), colnames(data))
  data <- matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames = dimnames) 
  
  G.R <- sapply(strsplit(colnames(data), "\\."), "[", 4 ) # "["(a, 4) almost a[4]
  
  
  if(is.null(k)||k == ""){ 
    ### you can choose the output data base on k
    return(NULL)
  } else {
    data.index <- grep(k, G.R)
    data_out <- data[, data.index]
    return(data_out) 
  }
}

mRNA_T_out <- splite_data(rt, "0..")
mRNA_N_out <- splite_data(rt, "1..")

mak_box_data <- function(mRNA_T_out, mRNA_N_out, gnam){
  ### prapare the data unit for boxplot
  g_tumor <- rbind("Tumor", mRNA_T_out[which(rownames(mRNA_T_out) == gnam), ])
  g_tumor <- data.frame(t(g_tumor), stringsAsFactors = F)
  g_normal <- rbind("Normal", mRNA_N_out[which(rownames(mRNA_N_out) == gnam), ])
  g_normal <- data.frame(t(g_normal), stringsAsFactors = F)
  data_unit <- rbind(g_tumor, g_normal)
  colnames(data_unit) <- c("group", gnam)
  data_unit[, 2] <- as.numeric(data_unit[, 2])
  return(data_unit)
}
data_unit <- mak_box_data(mRNA_T_out, mRNA_N_out, gnam)

mak_boxplot <- function(data_unit, gnam){
  T_Value <- t.test(data_unit[, 2]~group, data = data_unit, paired = FALSE)
  tiff(file = paste(cnam, gnam, ".tiff", sep = "" ), width = 512, height = 512)
  txt <- paste("p-value=", round(T_Value$p.value[1], digits = 4), sep = "")
  p = ggplot(data = data_unit, aes(data_unit[, 1], data_unit[, 2])) +
    xlab(cnam) +  ylab(paste("expression", "(log2FPKM)", sep = "")) +
    geom_boxplot(fill = c("blue", "red")) + ggtitle(txt) +
    theme_bw()+
    theme(plot.title = element_text(size=20, face="bold", hjust = 0.5), 
          axis.title.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          axis.text.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          axis.title.y = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
          axis.text.y = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
          panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
          panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5))
  print(p)
  dev.off()   
}
mak_boxplot(data_unit, gnam)

cnami <- "COAD"
gnami <- "hsa-miR-130b-5p"

miRNA_exp <- read.table(file = paste("/Users/stead/Desktop/GDC/", cnami, "/raw_data/", cnami, "_miRNA_HiSeq_gene", sep = ""), header = T, row.names = NULL, sep = "\t")
colnames(miRNA_exp) <- colnames(miRNA_exp)[-1]
miRNA_id_T <- read.table(file = "/Users/stead/Desktop/GDC/miRNA_id_transform.txt", header = T, sep = "\t")

splite_mi_data <- function(miRNA_exp, miRNA_id_T, k){
  G.R <- sapply(strsplit(colnames(miRNA_exp), "\\."), "[", 4 ) # "["(a, 4) almost a[4]
  if(is.null(k)||k == ""){ 
    ### TCGA-xx-xxxx-xx-xxxx. The sample id is the condition choosing the data_out
    return(NULL)
  } else {
    data.index <- grep(k, G.R) ## get the tumor group index
    miRNA_exp_T <- miRNA_exp[, data.index]
    miRNA_exp_T <- cbind(miRNA_exp[, 1], miRNA_exp_T)
    
    miRNA_id_T <- merge(miRNA_id_T, miRNA_exp_T, by = 1)
    miRNA_nam <- miRNA_id_T[, c(-1, -(3:6))]
    row.names(miRNA_nam) <- miRNA_nam[, 1]
    miRNA_T_out <- data.frame(apply(miRNA_nam[, -1], 2, as.numeric), stringsAsFactors = F)
    row.names(miRNA_T_out) <- miRNA_nam[, 1]
    return(miRNA_T_out)
  }
}



miRNA_T_out <- splite_mi_data(miRNA_exp, miRNA_id_T, "0.")


g_mRNA <- data.frame(mRNA_T_out[which(rownames(mRNA_T_out) == gnam), ])
g_mRNA <- cbind(row.names(g_mRNA), g_mRNA)
g_mRNA[, 1] <- gsub("(.*?)\\.(.*?)\\.(.*?)\\.(.*?)\\..*", "\\1\\.\\2\\.\\3\\.\\4", g_mRNA[, 1])
g_mRNA[, 1] <- substring(g_mRNA[, 1], 1, 15)

g_miRNA <- data.frame(t(miRNA_T_out[which(rownames(miRNA_T_out) == gnami), ]), stringsAsFactors = F)
g_miRNA <- data.frame(cbind(row.names(g_miRNA), g_miRNA))

data_unit <- merge(g_mRNA, g_miRNA, by = 1)
colnames(data_unit) <- c("sample_id", gnam, gnami)
data_unit$ENSG00000184863.9 <- as.numeric(data_unit$ENSG00000184863.9)
data_unit$`hsa-miR-130b-5p` <- as.numeric(data_unit$`hsa-miR-130b-5p`)

tiff(filename = paste(cnam, gnam, gnami, "correlation", ".tiff", sep = ""), width = 512, height = 512)
p = ggplot(data_unit, aes(data_unit[, 3], y =data_unit[, 2])) + 
  ylab(gnam) + xlab(gnami) +
  geom_point(size=3, colour = "grey60") + 
  stat_smooth(method=lm, level = 0.99) + theme_bw() +
  theme(plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.title.y = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 20, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 1),
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5))
print(p)
dev.off()  
