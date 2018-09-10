PathR <- "/Users/stead/Desktop/Lung_Cancer_CHQ_20170815/RNA_seq/20170901_data_analysis_batch1/QC/GCC"
setwd(PathR)
SamNam <- "abnormal.txt"
Fnam_Rt <- read.table(file = SamNam)
Fnam <- as.character(Fnam_Rt[, 1])

DrawData <- function(Fnam){
  Snam = unlist(strsplit(Fnam, ".", fixed = TRUE))[1]
  Rnam = unlist(strsplit(Fnam, ".", fixed = TRUE))[2]
  rt <- read.table(file = paste(Fnam, "_GC_content.txt", sep = ""), header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  rt_G <- rt[, c(1, 2)]
  #rt_G$Position <- factor(rt_G$Position, levels = unique(rt_G$Position))
  #rt_G$Base <- factor(rt_G$Base, levels = unique(rt_G$Base))
  rt_G$R <- rep(paste("R", Rnam, sep = ""), length(rt_G[, 1]))
  rt_G$Group <- rep(Snam, length(rt_G[,1]))
  return(rt_G)
}


###merge all samples data in one table
rt_M_list <- lapply(Fnam, DrawData)
rt_M <- data.frame(do.call(rbind, rt_M_list), stringsAsFactors = FALSE)

source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/theme_E.R")
tiff(filename = paste("GC_Content", "_", SamNam, ".tiff", sep = ""), width = 612, height = 512)
ggplot(rt_M, aes(x = GC, y = Content, group = Group, colour = Group)) + geom_line() + 
  ggtitle("GC content across all bases") + facet_grid(~R) +  theme_F +
  guides(colour = guide_legend(nrow = 2))
dev.off()
