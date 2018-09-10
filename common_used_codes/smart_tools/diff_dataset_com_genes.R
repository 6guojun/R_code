source('/Users/stead/Documents/SourceTree/R/common_used_codes/smart_tools/EA2SM.R')

rt_TCGA <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/', 'UVM', '/RNA_seq/TCGA-', 'UVM', '.htseq_fpkm.tsv', sep = ''), header = TRUE, 
                                 row.names = 1, sep = '\t', stringsAsFactors = FALSE)

rt_U133Plus2 <-  read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE37745', '/expression_data/', 'GSE37745', '_mat.txt', sep = '')
                                  , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

rt_U133A <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE14814', '/expression_data/', 'GSE14814', '_mat.txt', sep = '')
                       , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

rt_illu_WG6_3 <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE42127', '/expression_data/', 'GSE42127', '_mat.txt', sep = '')
                            , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_agilent_WHGM_4x44k <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE13213', '/expression_data/', 'GSE13213', '_mat.txt', sep = '')
                                    , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_agilent_UNC_4x44k <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE26939', '/expression_data/', 'GSE26939', '_mat.txt', sep = '')
                                   , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_affy_HG_U95 <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE83227', '/expression_data/', 'GSE83227', '_mat.txt', sep = '')
                             , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_agilent_HS_216K <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE11969', '/expression_data/', 'GSE11969', '_mat.txt', sep = '')
                                 , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_affy_H_FL <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', 'GSE68571', '/expression_data/', 'GSE68571', '_mat.txt', sep = '')
                           , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

###get symbol id
rt_TCGA_sym <-  EA2SM(rt_TCGA, in_type = 'ensembl_gene_id', out_type = 'hgnc_symbol', dup_type = 'mean')
rt_U133Plus2_sym <- EA2SM(rt_U133Plus2, in_type = 'affy_hg_u133_plus_2', out_type = 'hgnc_symbol', dup_type = 'mean')
rt_U133A_sym <- EA2SM(rt_U133A, in_type = 'affy_hg_u133a', out_type = 'hgnc_symbol', dup_type = 'mean')
rt_illu_WG6_3_sym <- EA2SM(rt_illu_WG6_3, in_type = 'illumina_humanwg_6_v3', out_type = 'hgnc_symbol', dup_type = 'mean')
#rt_agilent_WHGM_4x44k_sym <- EA2SM(rt_agilent_WHGM_4x44k, in_type = ' agilent_wholegenome_4x44k_v2', out_type = 'hgnc_symbol', dup_type = 'mean') ##**(id problem)
#rt_agilent_UNC_4x44k_sym <- EA2SM(rt_agilent_UNC_4x44k, in_type = ' agilent_wholegenome_4x44k_v2', out_type = 'hgnc_symbol', dup_type = 'mean') ##**(id problem)
#rt_affy_HG_U95_sym <- EA2SM(rt_affy_HG_U95, in_type = 'affy_hg_u95a', out_type = 'hgnc_symbol', dup_type = 'mean') ##**data set is null
#rt_agilent_HS_216K_sym <- EA2SM(rt_agilent_HS_216K, in_type = 'affy_hugenefl', out_type = 'hgnc_symbol', dup_type = 'mean') ##(id problem)
rt_affy_H_FL_sym <- EA2SM(rt_affy_H_FL, in_type = 'affy_hugenefl', out_type = 'hgnc_symbol', dup_type = 'mean')

###
setwd('/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/diff_platform_gene_symbol')
ensemble_U133Plus2_int <- data.frame(intersect(row.names(rt_TCGA_sym), row.names(rt_U133Plus2_sym)), stringsAsFactors = FALSE)
colnames(ensemble_U133Plus2_int) <- 'ensemble_U133Plus2_sym'
write.table(ensemble_U133Plus2_int, file = 'ensemble_U133Plus2_sym.txt', row.names = FALSE, col.names = TRUE, sep = '\t')

ensemble_U133A_int <- data.frame(intersect(row.names(rt_TCGA_sym), row.names(rt_U133A_sym)), stringsAsFactors = FALSE)
colnames(ensemble_U133A_int) <- 'ensemble_rt_U133A_sym'
write.table(ensemble_U133A_int, file = 'ensemble_U133A_sym.txt', row.names = FALSE, col.names = TRUE, sep = '\t')

ensemble__U133Plus2_U133A_int <- data.frame(intersect(row.names(rt_U133Plus2_sym), intersect(row.names(rt_TCGA_sym), row.names(rt_U133A_sym))), stringsAsFactors = FALSE)
colnames(ensemble__U133Plus2_U133A_int) <- 'ensemble__U133Plus2_U133A_sym'
write.table(ensemble__U133Plus2_U133A_int, file = 'ensemble__U133Plus2_U133A_sym.txt', row.names = FALSE, col.names = TRUE, sep = '\t')

###add other platform
ensemble__U133Plus2_U133A_illu_WG6_3_int <- data.frame(intersect(row.names(rt_illu_WG6_3_sym),intersect(row.names(rt_U133Plus2_sym), intersect(row.names(rt_TCGA_sym), row.names(rt_U133A_sym)))), 
                                            stringsAsFactors = FALSE)
colnames(ensemble__U133Plus2_U133A_illu_WG6_3_int) <- 'ensemble__U133Plus2_U133A_illu_WG6_3_sym'
write.table(ensemble__U133Plus2_U133A_illu_WG6_3_int, file = 'ensemble__U133Plus2_U133A_illu_WG6_3_sym.txt', row.names = FALSE, col.names = TRUE, sep = '\t')

