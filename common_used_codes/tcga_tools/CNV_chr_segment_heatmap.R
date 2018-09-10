setwd('/Users/stead/Desktop/KIRP_submit/test')
rt_uvm_cnv <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/copy_num/TCGA-", cancer , ".masked_cnv.tsv", sep = ""),
                         sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rt_uvm_cnv$arm <- NA
rt_uvm_cnv$n.probes <- NA
rt_uvm_cnv_c <- rt_uvm_cnv[, c(1, 2, 6, 3, 4, 7, 5)]
colnames(rt_uvm_cnv_c) <- c("sampleID", "chrom", "arm", "start.pos", "end.pos", 'n.probes', "mean")
sam_cnv <- matrix(unlist(strsplit(rt_uvm_cnv_c$sampleID, '-')), ncol = 4, byrow = TRUE)
rt_uvm_cnv_c$sampleID <- paste(sam_cnv[, 1], sam_cnv[, 2], sam_cnv[, 3], sep = "-")

###get cnv and exp common samples id
int_cnv_id <- intersect(rt_uvm_cnv_c$sampleID, row.names(rt_mvm_sur_all_order))
rt_mvm_sur_cnv <- rt_mvm_sur_all_order[match(int_cnv_id, row.names(rt_mvm_sur_all_order), nomatch = 0), ]
rt_mvm_sur_cnv_order <- rt_mvm_sur_cnv[order(rt_mvm_sur_cnv$risk_score), ]
pos_cnv <- match(rt_uvm_cnv_c$sampleID, row.names(rt_mvm_sur_cnv_order), nomatch = 0)
cnv_order <- order(pos_cnv, decreasing = TRUE)
rt_uvm_cnv_o <- rt_uvm_cnv_c[cnv_order, ]

pdf(file = 'cnv_risk.pdf', 8, 10)
plotHeatmap(segments = rt_uvm_cnv_o, upper.lim = 1, lower.lim = -1, chrom = NULL, colors = c('purple3', 'white', 'red'))
dev.off()