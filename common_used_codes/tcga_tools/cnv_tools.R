###CNV heatmap

rt_uvm_cnv <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/UVM/copy_num/TCGA-UVM.masked_cnv.tsv", sep = '\t', header = TRUE)
rt_uvm_cnv$arm <- NA
rt_uvm_cnv$n.probes <- NA
rt_uvm_cnv_c <- rt_uvm_cnv[, c(1, 2, 6, 3, 4, 7, 5)]
colnames(rt_uvm_cnv_c) <- c("sampleID", "chrom", "arm", "start.pos", "end.pos", 'n.probes', "mean")
plotHeatmap(segments = rt_uvm_cnv_c, upper.lim = 1, lower.lim = -1, chrom = NULL, colors = c('purple3', 'white', 'red'))


