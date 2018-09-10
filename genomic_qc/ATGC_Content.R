library(dplyr)
library(tidyr)
library(ggplot2)
PathR <- "/Users/stead/Desktop/Lung_Cancer_CHQ_20170815/RNA_seq/20170901_data_analysis_batch1/QC/ATGC"
setwd(PathR)
SamNam <- "abnormal.txt"
Fnam_Rt <- read.table(file = SamNam)
Fnam <- as.character(Fnam_Rt[, 1])

###Sort data
DrawData <- function(Fnam){
  Snam = unlist(strsplit(Fnam, ".", fixed = TRUE))[1]
  Rnam = unlist(strsplit(Fnam, ".", fixed = TRUE))[2]
  rt <- read.table(file = paste(Fnam, "_base_sequence_content.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  rt_G <- gather(rt, Base, Value, -Base)
  colnames(rt_G) <- c("Position", "Base", "Value")
  rt_G$Position <- factor(rt_G$Position, levels = unique(rt_G$Position))
  rt_G$Base <- factor(rt_G$Base, levels = unique(rt_G$Base))
  rt_G$R <- rep(paste("R", Rnam, sep = ""), length(rt_G[, 1]))
  rt_G$Group <- rep(Snam, length(rt_G[,1]))
  return(rt_G)
}
   
###merge all samples data in one table
rt_M_list <- lapply(Fnam, DrawData)
rt_M <- data.frame(do.call(rbind, rt_M_list), stringsAsFactors = FALSE)

source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/Theme_F.R")
tiff(filename = paste("ATCG_Content", SamNam, ".tiff", sep = ""), width = 512, height = 512)
ggplot(rt_M, aes(x = Position, y = Value, group = Group, colour = Group), size = 5) + geom_line() + 
  facet_grid(R~Base, space="free_x", scales="free_x") + 
   xlab("Position in read (bp)") + ylab("Percent (%)") +
  ggtitle("Sequence content  across all bases") + theme_F +
  guides(colour = guide_legend(nrow = 2))
dev.off()

