library("ggplot2")
rt <- read.table(file = "/Users/stead/Desktop/GDC/LUAD/merge_result/27_log2exp_pac_heatmap/cancer_T_genes_merge.txt", header = T, sep = "\t", stringsAsFactors = F)
a <- c("ENSG00000019169.10", "ENSG00000122852.13", "ENSG00000168484.11", "ENSG00000185303.14", 
       "ENSG00000203878.10", "ENSG00000231322.4", "ENSG00000260695.1", "ENSG00000248608.2",
       "ENSG00000168878.15", "ENSG00000203878.10", "ENSG00000112175.7", "ENSG00000166961.13")
b <- c("MARCO", "SFTPA1", "SFTPC", "SFTPA2", "CHIAP2", "RPL13AP17", "RP11-513N24.11", "RP11-206P5.2", 
       "SFTPB", "CHIAP2", "BMP5", "MS4A15")

#d <- c("green1", "orangered1", "yellowgreen", "cyan1", "cyan2", "pink1", "orchid4", 
#       "orchid3", "orchid2", "orchid1", "pink2", "steelblue2",  "steelblue1", 
#       "brown", "mediumpurple2", "mediumpurple3","mediumpurple4",  "hotpink1", "god3", "hotpink3", 
#       "skyblue1", "skyblue2", "skyblue3", "green1", "tan1", "tan2", "tan3")

mak_merge <- function(k, rt){
  rt_gnam <- cbind(b[k], rt[, "group"], rt[, a[k]])
  return(rt_gnam)
}
k <- c(1:length(a))
rt_out <- lapply(k, mak_merge, rt)
### you can choose which genes you want to add 
rt_out_d <- do.call("rbind", rt_out[7:12])
rt_box <- data.frame(rt_out_d, stringsAsFactors = F)
colnames(rt_box) <- c("gnam", "group", "exp")
#color_G <- factor(rt_box$group, labels = d)


rt_box$exp <- as.numeric(rt_box$exp)
rt_box$gnam <- factor(rt_box$gnam, levels = b[7:12])

filename <- "SUR_T_boxplot"
mak_facet_box <- function(rt_box){
  tiff(file = paste(filename, ".tiff", sep = ""), width = 4000, height = 1800)
  p = ggplot(data = rt_box, aes(group,  rt_box$exp)) +
   facet_wrap(~gnam, ncol = 3, nrow = 2, drop = T) + ###if you gnam is a factor, the facet will order of the factor
  # facet_grid(~ gnam)+
    xlab("cancer_type") +  ylab("gene expression (log2FPKM)")+
    geom_boxplot(aes(fill= group)) +
    theme_bw()+ 
    theme(plot.title = element_text(size=60, face="bold", hjust = 0.5), 
          strip.text = element_text(size = 60),
          axis.title.x = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 1),
          axis.text.x = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 90),
          axis.title.y = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 450),
          axis.text.y = element_text(size = 50, color = "black",  vjust = 0.5, hjust = 0.5, angle = 1),
          panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
          panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5),
          legend.title = element_text(size=40, color="black"),
          legend.key.size = unit(2, "cm"),
          legend.text = element_text(colour = 'black', size = 30))
  print(p)
  dev.off() 
}

mak_facet_box(rt_box)
  
 