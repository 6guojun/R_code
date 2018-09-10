library(ggplot2)
setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/mapping")
rt <- read.table(file = "reads_mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(rt) <- c("samples", "total_reads", "mapping_rate")


rt_sam <- read.table(file = "samples.txt", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)

rt$group <- rt_sam$group
rt$samples <- rt_sam$samples
rt$color <- rt_sam$color
rt_reads <- rt[, c("samples", "total_reads", "group", "color")]
rt_mapping <- rt[, c("samples", "mapping_rate", "group", "color")]

rt_reads$group <- factor(rt_reads$group, levels = unique(rt_reads$group), labels = unique(rt_reads$group))
rt_reads <- rt_reads[order(rt_reads$total_reads), ]
rt_reads$total_reads <- round(rt_reads$total_reads, digits = 3)
rt_reads$samples <- factor(rt_reads$samples, levels = rt_reads$samples)
source("/Users/stead/Documents/SourceTree/R/ggplot/ggplot_theme/theme_E.R")
tiff(filename = "total_reads.tiff", width = 1024, height = 512)
ggplot(rt_reads, aes(x = samples, y = total_reads, fill = group)) + 
  geom_bar(stat="identity") + ggtitle("total_reads") + scale_fill_manual(values = unique(rt$color)) +
  theme_E
dev.off()

rt_mapping$group <- factor(rt_mapping$group, levels = unique(rt_mapping$group), labels = unique(rt_mapping$group))
rt_mapping <- rt_mapping[order(rt_mapping$mapping_rate), ]
rt_mapping$samples <- factor(rt_mapping$samples, levels = rt_mapping$samples)
tiff(filename = "mapping_reate.tiff", width = 1024, height = 512)
ggplot(rt_mapping, aes(samples, mapping_rate, fill = group)) + scale_fill_manual(values = unique(rt$color)) +
  geom_bar(stat="identity") + ylab("mapping_rate (%)") + ggtitle("mapping_rate") +
  theme_E
dev.off()
 
  

