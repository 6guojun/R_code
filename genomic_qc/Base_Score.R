PathR <- "/Users/stead/Desktop/Lung_Cancer_CHQ_20170815/RNA_seq/20170901_data_analysis_batch1/QC/BSQ"
setwd(PathR)
SamNam <- "xaf"
Fnam_Rt <- read.table(file = SamNam)
Fnam <- as.character(Fnam_Rt[, 1])


###Sort data
DrawData <- function(Fnam){
  Snam = unlist(strsplit(Fnam, ".", fixed = TRUE))[1]
  Rnam = unlist(strsplit(Fnam, ".", fixed = TRUE))[2]
  rt <- read.table(file = paste(Fnam, "_base_sequence_quality.txt", sep = ""), header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  rt_G <- rt[, c("Base", "Mean")]
  rt_G$Base <- factor(rt_G$Base, levels = unique(rt_G$Base))
  rt_G$R <- rep(paste("R", Rnam, sep = ""), length(rt_G[, 1]))
  rt_G$Group <- rep(Snam, length(rt_G[,1]))
  return(rt_G)
}

###merge all samples data in one table
rt_M_list <- lapply(Fnam, DrawData)
rt_M <- data.frame(do.call(rbind, rt_M_list), stringsAsFactors = FALSE)

source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/theme_E.R")
tiff(filename = paste("Base_score", "-", SamNam, ".tiff", sep = ""), width = 720, height = 512)
ggplot(rt_M, aes(x = Base, y = Mean,  group = Group, colour = Group)) + geom_line() + 
  facet_grid(~R) + xlab("Position in read (bp)") + ylab("Scores (mean)") +
  ggtitle("Quality scores across all bases") +  theme_E
dev.off()

