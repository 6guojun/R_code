library(GEOquery)

gse81089 <- getGEO('GSE81089', GSEMatrix = TRUE)
rt_cli_e81089 <- pData(phenoData(gse81089[[1]]))
rt_cli_e81089_S <- rt_cli_e81089[, c("age:ch1", "dead:ch1", "gender:ch1", "stage tnm:ch1", "surgery date:ch1",  "histology:ch1", "tumor (t) or normal (n):ch1", "vital date:ch1")]
colnames(rt_cli_e81089_S) <- c('age', 'OS_status', 'gender', 'stage', 'surgery_date', 'histology', 'tissue_type', 'vital_data')

rt_cli_e81089_S_AD <- rt_cli_e81089_S[which(rt_cli_e81089_S$histology == "2"), ]
write.table(rt_cli_e81089_S, file = 'raw_clinical_data.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
rt_sur <- read.table(file = '/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/GSE81089/survival_data/survival_data.txt', row.names = 1, 
           header = TRUE, stringsAsFactors = FALSE)
rt_sur_AC <- rt_sur[match(row.names(rt_cli_e81089_S_AD), row.names(rt_sur), nomatch = 0), ]
write.table(rt_sur_AC, file = 'survival_data.txt', row.names = TRUE, col.names = TRUE, sep = '\t')

rt_exp <- read.table(file = '/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/GSE81089/expression_data/GSE81089_FPKM_cufflinks.tsv', header = TRUE, row.names = 1, 
                     sep = '\t', stringsAsFactors = FALSE)
colnames(rt_exp) <- row.names(rt_cli_e81089_S)[match(colnames(rt_exp), rt_cli_e81089_S$tissue_type, nomatch = 0)]

rt_exp <- data.frame(apply(rt_exp, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
write.table(rt_exp, file = 'GSE81089_mat.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
rm(rt_exp) 

GSEID <- 'GSE50081'
setwd(paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', GSEID, '/expression_data', sep = ""))
rt_exp <- read.table(file = paste(GSEID, '_mat.txt', sep = ""), header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
rt_exp_log <- data.frame(apply(rt_exp, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
write.table(rt_exp_log, file = paste(GSEID, '_mat_log.txt', sep = ""), row.names = TRUE, col.names = TRUE, sep = '\t')
rm(rt_exp, rt_exp_log)

GSEID <- 'GSE3141'
setwd(paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', GSEID, '/expression_data', sep = ""))
rt_exp <- read.table(file = paste(GSEID, '_mat.txt', sep = ""), header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
rt_exp_log <- data.frame(apply(rt_exp, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
write.table(rt_exp_log, file = paste(GSEID, '_mat_log.txt', sep = ""), row.names = TRUE, col.names = TRUE, sep = '\t')
rm(rt_exp, rt_exp_log)

GSEID <- 'GSE31210'
setwd(paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', GSEID, '/expression_data', sep = ""))
rt_exp <- read.table(file = paste(GSEID, '_mat.txt', sep = ""), header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
rt_exp_log <- data.frame(apply(rt_exp, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
write.table(rt_exp_log, file = paste(GSEID, '_mat_log.txt', sep = ""), row.names = TRUE, col.names = TRUE, sep = '\t')
rm(rt_exp, rt_exp_log)

GSEID <- 'GSE19188'
setwd(paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', GSEID, '/expression_data', sep = ""))
rt_exp <- read.table(file = paste(GSEID, '_mat.txt', sep = ""), header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
rt_exp_log <- data.frame(apply(rt_exp, 2, function(x){log2(x + 1)}), stringsAsFactors = FALSE)
write.table(rt_exp_log, file = paste(GSEID, '_mat_log.txt', sep = ""), row.names = TRUE, col.names = TRUE, sep = '\t')
rm(rt_exp, rt_exp_log)


gse22138 <- getGEO('GSE22138', GSEMatrix = TRUE)
rt_cli_e22138 <- pData(phenoData(gse22138[[1]]))

gse13213 <- getGEO('GSE13213', GSEMatrix = TRUE)
rt_cli_e13213 <- pData(phenoData(gse13213[[1]]))

GSE83227 <- getGEO('GSE83227', GSEMatrix = TRUE)
GSE83227_cli <- pData(phenoData(GSE83227[[1]]))

GSE10245 <- getGEO('GSE10245', GSEMatrix = TRUE)
GSE10245_cli <- pData(phenoData(GSE10245[[1]]))

GSE11969 <- getGEO('GSE11969', GSEMatrix = TRUE)
GSE11969_cli <- pData(phenoData(GSE11969[[1]]))

GSE10445 <- getGEO('GSE10445', GSEMatrix = TRUE)
GSE10445_cli <- pData(phenoData(GSE10445[[1]]))

GSE18842 <- getGEO('GSE18842', GSEMatrix = TRUE)
GSE18842_cli <- pData(phenoData(GSE18842[[1]]))

GSE12667 <- getGEO('GSE12667', GSEMatrix = TRUE)
GSE12667_cli <- pData(phenoData(GSE18842[[1]]))

GSE33356 <- getGEO('GSE33356', GSEMatrix = TRUE)
GSE33356_cli <- pData(phenoData(GSE18842[[1]]))


GSE28571 <- getGEO('GSE28571', GSEMatrix = TRUE)
GSE28571_cli <- pData(phenoData(GSE28571[[1]]))

GSE77803 <- getGEO('GSE77803', GSEMatrix = TRUE)
GSE77803_cli <- pData(phenoData(GSE77803[[1]]))

GSE63074 <- getGEO('GSE63074', GSEMatrix = TRUE)
GSE63074_cli <- pData(phenoData(GSE63074[[1]]))

GSE5843 <- getGEO('GSE5843', GSEMatrix = TRUE)
GSE5843_cli <- pData(phenoData(GSE5843[[1]]))

GSE7339 <- getGEO('GSE7339', GSEMatrix = TRUE)
GSE7339_cli <- pData(phenoData(GSE7339[[1]]))

GSE120622 <- getGEO('GSE120622', GSEMatrix = TRUE)
GSE120622_cli <- pData(phenoData(GSE120622[[1]]))
GSE120622_cli_2 <- pData(phenoData(GSE120622[[2]]))
write.table(GSE120622_cli_2, file = 'clinical_data_raw.txt', row.names = TRUE, col.names = TRUE, sep = '\t')

GSE10245 <- getGEO('GSE10245', GSEMatrix = TRUE)
GSE10245_cli <- pData(phenoData(GSE10245[[1]]))
GSE10245_cli_2 <- pData(phenoData(GSE10245[[2]]))

GSE1037 <- getGEO('GSE1037', GSEMatrix = TRUE)
GSE1037_cli <- pData(phenoData(GSE1037[[1]]))
GSE1037_cli_2 <- pData(phenoData(GSE1037[[2]]))

GSE5123 <- getGEO('GSE5123', GSEMatrix = TRUE)
GSE5123_cli <- pData(phenoData(GSE5123[[1]]))
GSE5123_cli_2 <- pData(phenoData(GSE5123[[2]]))

GSE4573 <- getGEO('GSE4573', GSEMatrix = TRUE)
GSE4573_cli <- pData(phenoData(GSE4573[[1]]))
GSE4573_cli_2 <- pData(phenoData(GSE4573[[2]]))
