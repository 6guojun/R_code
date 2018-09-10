library(ggplot2)
rt <- read.table(file = "ARD_LCL6_QS.txt", header = TRUE, sep = "\t")
tiff(filename = "Genomic_QS_1.tiff", width = 1024, height = 712)
ggplot(rt, aes(x = CYCLE, y = MEAN_QUALITY), ) + geom_line(colour = "blue", size = 1.5) + 
   ylim(20, 50) + theme_A
dev.off()

