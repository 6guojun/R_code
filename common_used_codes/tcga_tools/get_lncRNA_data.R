setwd('/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/nocoding_RNA')
rt_lnc <- read.table(file = '/Users/stead/Downloads/gencode.v28.long_noncoding_RNAs.gtf', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
rt_lnc_gene <- rt_lnc[which(rt_lnc$V3 == 'gene'), ]
lnc_list <- strsplit(rt_lnc_gene[, 9], ";")
lnc_mat <- do.call(rbind, lnc_list)
lnc_id <- strsplit(lnc_mat[, 1], ' ')
lnc_type <- strsplit(lnc_mat[, 2], ' ')
lnc_name <- strsplit(lnc_mat[, 3], ' ')

lnc_dat <- cbind(do.call(rbind, lnc_id)[, 2], do.call(rbind, lnc_name)[, 3], do.call(rbind, lnc_type)[, 3])
colnames(lnc_dat) <- c('ens_id', 'gene_name', 'lnc_type')
write.table(lnc_dat, file = 'lncRNA_ensid_name_type.txt', row.names = TRUE, col.names = TRUE, sep = '\t')

lnc <- read.table(file = '/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/lncRNA_ensid_name_type.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
cancer1 <- 'UVM'
rt <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer1, "/RNA_seq/TCGA-", cancer1, ".htseq_fpkm.tsv", sep = ""), 
                    sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt_lnc <- merge(lnc, rt, by = 1)

write.table(rt_lnc, file = paste(cancer1, '_lncRNA_data.txt'), sep = '\t', row.names = FALSE, col.names = TRUE)

###pseudogene
rt_pse <- read.table(file = '/Users/stead/Downloads/gencode.v7.pseudogene.txt', header = FALSE, sep = '\t', stringsAsFactors = FALSE)
pse_list <- strsplit(rt_pse[, 9], ";")
pse_mat <- do.call(rbind, pse_list)
pse_id <- strsplit(pse_mat[, 1], ' ')
pse_type <- strsplit(pse_mat[, 3], ' ')
pse_name <- strsplit(pse_mat[, 5], ' ')

pse_dat <- cbind(do.call(rbind, pse_id)[, 2], do.call(rbind, pse_name)[, 3], do.call(rbind, pse_type)[, 3])
colnames(pse_dat) <- c('ens_id', 'gene_name', 'pse_type')
write.table(pse_dat, file = 'pseudogene_ensid_name_type.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
