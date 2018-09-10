library(ggplot2)
library(devtools)
library(easyGgplot2)
library(tidyr)
library(dplyr)

setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_batch4/fastq_screen")
SamNam <- "samples.txt"
Fnam_Rt <- read.table(file = SamNam)
Fnam <- as.character(Fnam_Rt[, 1])

###sort the data
DrawData <- function(Fnam){
  rt <- read.table(file = paste(Fnam, "_FQScreen.txt", sep = ""), header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  rt$Hit_no_genomes <- c(rep(0, 11), rt[12, 2])
  rt <- rt[, -2]
  rt[12, c(2:5)] <- c(rep(0, 4))
  rt_G <- gather(rt, group, mapped, -Genome)
  rt_G$Snam <- rep(Fnam, length(rt_G[, 1]))
  return(rt_G)
}

rt_M_list <- lapply(Fnam, DrawData)
rt_M <- data.frame(do.call(rbind, rt_M_list), stringsAsFactors = FALSE)
rt_M$Genome <- factor(rt_M$Genome, levels = rev(unique(rt_M$Genome)))

source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/theme_F.R")
rt_M$Snam <- factor(rt_M$Snam, levels = as.character(unique(rt_M$Snam)))
tiff(filename = paste("FQS", "-", SamNam, ".tiff", sep = ""), width = 880, height = 512)
ggplot(rt_M, aes(fill = group, y = mapped, x = Genome)) +  
  facet_grid(~Snam) + geom_bar( stat="identity") + theme_F + 
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
dev.off()

