library('biomaRt')

setwd('/Users/stead/Desktop/subtype_analysis/signature/LUAD/arrary')
GSE_ID = 'GSE37745'
coef_value <- read.table(file = '/Users/stead/Desktop/subtype_analysis/signature/LUAD/LUAD_0.6/LUAD_genes_coef.txt', header = TRUE, row.names = 1, sep = '\t', 
                         stringsAsFactors = FALSE)
uvm_risk_genes <- read.table(file = '/Users/stead/Desktop/subtype_analysis/signature/LUAD/all/uvm_elements_uvm_cox.txt', header = TRUE, row.names = 1, sep = '\t', 
                             stringsAsFactors = FALSE)

##arry data
rt_exp <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', GSE_ID, '/expression_data/', GSE_ID, '_mat.txt', sep = '')
                     , sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rt_sur <- read.table(file = paste('/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/LUAD/', GSE_ID, '/survival_data/survival_data.txt', sep = ''), 
                     header = TRUE, sep = '\t', stringsAsFactors = FALSE)

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
eg_arry <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol", "ensembl_gene_id"), values = row.names(rt_exp), mart = ensembl)

eg_arry_d <- eg_arry[-which(eg_arry$affy_hg_u133_plus_2 == ""), ]


ens_pos <- match(coef_value$risk_gens, eg_arry$ensembl_gene_id, nomatch = 0)
sym_pos <- match(row.names(uvm_risk_genes)[1:11], eg_arry$hgnc_symbol, nomatch = 0)
arry_pos <- match(eg_arry$affy_hg_u133_plus_2[ens_pos], row.names(rt_exp), nomatch = 0)

###count risk score
rt_exp_risk <- rt_exp[arry_pos, ]
rt_risk_score_mat <- rt_exp_risk*(coef_value$coef_value[-which(eg_arry$affy_hg_u133_plus_2[ens_pos] == "")])
risk_score_mat <- data.frame(cbind(colnames(rt_risk_score_mat), apply(rt_risk_score_mat, 2, sum)), stringsAsFactors = FALSE)
colnames(risk_score_mat) <- c('samples_id', 'risk_score')

###prepare survival anaysis data
exp_sur_int <- intersect(risk_score_mat$samples_id, rt_sur$Accession)
rt_sur_m <- rt_sur[match(exp_sur_int, rt_sur$Accession, nomatch = 0), ]
rt_score_mat_m <- risk_score_mat[match(exp_sur_int, risk_score_mat$samples_id, nomatch = 0), ]
rt_sur_risk <- data.frame(cbind(rt_sur_m, rt_score_mat_m$risk_score), stringsAsFactors = FALSE)
colnames(rt_sur_risk) <- c('samples_id', 'OS_Status', 'OS_Time', 'risk_score')
rt_sur_risk$OS_Status <- gsub('ALIVE', 0, rt_sur_risk$OS_Status)
rt_sur_risk$OS_Status <- gsub('DEAD', 1, rt_sur_risk$OS_Status)
rt_sur_risk$risk_score <- as.numeric(as.character(rt_sur_risk$risk_score))

#do survival analysis
source('/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/survival_analysis_2.R')
source('/Users/stead/Documents/SourceTree/R/common_used_codes/survival_analysis/survival_miner.R')
source('/Users/stead/Documents/SourceTree/R/common_used_codes/survival_analysis/ggsurv_theme/ggsurv_theme_A.R')
Rt_Exp_Cli <- rt_sur_risk
Gnam <- 'risk_score'
Rt_Exp_Cli[, Gnam][Rt_Exp_Cli[, Gnam] > median(Rt_Exp_Cli[, Gnam])] <- 'high_risk'
Rt_Exp_Cli[, Gnam][which(Rt_Exp_Cli[, Gnam] != 'high_risk')] <- 'low_risk'
Rt_Exp_Cli$OS_Time <- as.numeric(as.character(Rt_Exp_Cli$OS_Time))
Rt_Exp_Cli$OS_Status <- as.numeric(as.character(Rt_Exp_Cli$OS_Status))
fit <- survfit(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli[, Gnam], data = Rt_Exp_Cli)
ggsurvplot(fit, pval = TRUE)

