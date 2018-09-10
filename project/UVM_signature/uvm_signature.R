library(plyr)
library(ggplot2)
library(glmnet)
library(survival)
library(survminer)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(maftools)
library(survivalROC)
library(copynumber)

cancer <- "UVM"

setwd(paste("/Users/stead/Desktop/subtype_analysis/signature/", cancer, sep = ''))
source("/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/get_exp_sur_cli_merge_data.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/get_tcga_sample_type.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/clinical_data_sort.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/cox_analysis.R")


rt_exp <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/RNA_seq/TCGA-", cancer, ".htseq_fpkm.tsv", sep = ""),
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rt_cli <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/phenotype/TCGA-", cancer, ".GDC_phenotype.tsv", sep = ""),
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE, fill = TRUE, quote = "", na.strings = "NA")
rt_sur <- read.table(file = paste("/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/", cancer, "/phenotype/TCGA-", cancer, ".survival.tsv", sep = ""),
                     sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)


rt_T <- split_tcga_tn(rt_exp[, -1], sam_type = "tumor")


###get merged exp, cli and sur
cancer_list <- GetExpSurCli(rt_T, rt_cli, rt_sur)
rt_T_m <- cancer_list[[1]]
rt_cli_m <- cancer_list[[2]]
rt_sur_m <- cancer_list[[3]]
cancer_cli_sort <- CliSort(rt_cli_m, cancer)
rt_cli_m <- cancer_cli_sort

###get cox table which can be used to perform cox analysis
exp_d_pos <- apply(rt_T_m, 1, function(x){median(x) > 0})
rt_T_m_d <-rt_T_m[exp_d_pos, ]
rt_T_m_d <- apply(rt_T_m_d, 2, function(x){log2(x + 1)})
unam <- rt_exp$Ensembl_ID[exp_d_pos] 

rt_cancer <- data.frame(cbind(t(rt_T_m_d), rt_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)
colnames(rt_cancer) <- c(unam, "OS_Status", "OS_Time")

###one element analysis
uvm_list <- lapply(unam, uvm_count, rt_cancer)
data_uvm <- do.call(rbind, uvm_list)

###mutiple cox analysis
sig_genes <-  row.names(data_uvm)[data_uvm$pvalue < 0.01]
rt_mvm <- rt_T_m[match(sig_genes, rt_exp$Ensembl_ID, nomatch = 0), ]
row.names(rt_mvm) <- sig_genes
rt_mvm_sur <- data.frame(cbind(t(rt_mvm), rt_sur_m[, c("X_OS_IND", "X_OS")]), stringsAsFactors = FALSE)

###glmnet
source('/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/glamnet_cox_modle.R')
x <- data.matrix(t(rt_mvm))
y <- data.matrix(rt_mvm_sur[, c('X_OS', 'X_OS_IND')])
colnames(y) <- c("time", "status")
x_coef_value_list <- GLCoxMain(2000, x, y, 10, 6, theme_E)

risk_genes <- strsplit(x_coef_value_list[1][[1]][[1]], ';')[[1]]
x_coef<- x_coef_value_list[2][[1]][[1]]
coef_value_d <- x_coef_value_list[3][[1]][[1]]
GC_mat <- x_coef_value_list[4][[1]]

###get uvm risk genes table 
data_uvm_risk_gens <- data_uvm[risk_genes, ]
write.table(data_uvm_risk_gens, file = 'uvm_risk_gens.txt', row.names = TRUE, col.names = TRUE)

#do frequency boxplot
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_E.R")
pdf(file = paste("fren_boxplot2.pdf", sep = ""), width = 10, height = 6)
p <- ggplot(GC_mat, aes(x = reorder(type, freuency) , y = freuency)) + 
  geom_bar(stat = "identity", fill = 'cyan3')  + labs(title = 'the frequency of models') + ylab("frequency")  + xlab('type') +
  geom_text(aes(label = freuency),vjust = 0.5, hjust = 0.5, color = 'black', size = 5) + theme_E
print(p)
dev.off()



###make a cox risk modle

MakSum <- function(i){
  coef_exp <- x_coef[, 1]*coef_value_d[1]
  for(i in 2 : i){
    coef_expi <- x_coef[, i]*coef_value_d[i]
    coef_exp <- coef_exp + coef_expi
  }
  return(coef_exp)
}

i = length(coef_value_d)
coef_score <- MakSum(i)


#survival analysis
risk_genes_ens <- matrix(unlist(strsplit(risk_genes, '[.]')), ncol = 2, byrow = TRUE)[, 1]
risk_genes_sym <- bitr(risk_genes_ens, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db', drop = FALSE)[, 'SYMBOL']
risk_genes_sym <- c('AL137784.1', "S100A13", "CA12", "MGLL", "PARP8", "FAM189A2", "ZBED1", "MIR4655", "ZNF497",  'AC023790.2', "RNF208", "TCTN1", 'FABP5P1', "GRIN2A", 'AC092821.1', "SIRT3", "MMP9", 'AC010442.3')
rt_mvm_sur_min <- rt_mvm_sur[, c(risk_genes, 'X_OS', 'X_OS_IND')]
colnames(rt_mvm_sur_min) <- c(risk_genes_sym, 'X_OS', 'X_OS_IND')

rt_mvm_sur_min$score <- coef_score
colnames(rt_mvm_sur_min)[(length(colnames(rt_mvm_sur_min)) -2) : length(colnames(rt_mvm_sur_min))] <- c("OS_Time", "OS_Status", 'risk_score')


source('/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/survival_analysis_2.R')
gnam <- colnames(rt_mvm_sur_min)[-match(c("OS_Time", "OS_Status"), colnames(rt_mvm_sur_min))]
lapply(gnam, SurMainFun, rt_mvm_sur_min, DatType = "ConType", feature = "OS", OutType = "SurF", SurvType = 'ALL', len_a = 'high_exp', len_b = 'low_exp')


###time dependent ROC 
source('/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/time_ROC.R')
CountTimeROC(rt_mvm_sur_min, 'OS_Time', 'OS_Status', 'risk_score', predict_time = 1825, lambda = 0.05, 5)  


###do clinical and genomic analysis
rt_mvm_sur_min$stage <- rt_cli_m$stage
rt_mvm_sur_min$gender <- rt_cli_m$gender
rt_mvm_sur_min$recurrence <- rt_cli_m$recurrence
rt_mvm_sur_min$age <- rt_cli_m$age

###do complex heatmap
color_type_table <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/color.txt", sep = "\t", stringsAsFactors = FALSE)
color_type <- color_type_table$V2

order_score <- order(rt_mvm_sur_min$risk_score)
rt_mvm_sur_order <- rt_mvm_sur_min[order_score, ] 
os_time <- rt_mvm_sur_order$OS_Time
status <- as.character(rt_mvm_sur_order$OS_Status)
risk_score <- rt_mvm_sur_order$risk_score
stage <- rt_mvm_sur_order$stage
gender <- rt_mvm_sur_order$gender
recurrence <- rt_mvm_sur_order$recurrence
age <- rt_mvm_sur_order$age
exp <- apply(rt_mvm_sur_order[, risk_genes_sym], 1, scale)
row.names(exp) <- risk_genes_sym


ha1 = HeatmapAnnotation(recurrence = recurrence, gender = gender, stage = stage, 
                        col = list(recurrence= structure(names = unique(recurrence), c('grey', 'black', 'white')),
                        gender = structure(names = unique(gender), c('grey', 'black')),
                        stage = structure(names = unique(stage), c('black', 'grey', 'white'))), 
                        risk_score = anno_points(as.numeric(risk_score), gp = gpar(col = ifelse( risk_score > median(risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE), 
                        os_time =  anno_points(as.numeric(os_time), gp = gpar(col = ifelse( status == 1, "black", "grey")),width = unit(4, "cm"), axis = TRUE),
                        show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                        annotation_height = unit(c(5, 5, 5, 30, 30), "mm"), show_annotation_name = TRUE)


pdf(file = paste(cancer, "_MCP_subtype", "_cluster.pdf", sep = ""), 10, 8)
ht_list <- Heatmap(exp, name = "center scaled expression", col = colorRamp2(c(min(exp) , 0, max(exp)), c(" blue", "white", "red")),
                  bottom_annotation = ha1, bottom_annotation_height = unit(10, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = FALSE,
                   show_row_names = TRUE, show_column_names = FALSE) 
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()


###mvm cox analysis
rt_mvm_cox <- rt_mvm_sur_order[-c(grep(3, rt_mvm_sur_order$stage), grep(3,rt_mvm_sur_order$recurrence)), ]
write.table(rt_mvm_cox, file = 'risk_mvm_cox.txt', row.names = TRUE, col.names = TRUE, sep = '\t')
con_nam <- "risk_score"
log_nam <- c("stage", "gender", "recurrence")
data_mvm_cox  <- mvm_count(con_nam, log_nam, 'UVM_MVM', rt_mvm_cox)

#any stage survival analysis
Rt_Exp_Cli_stageI <- rt_mvm_sur_order[which(rt_mvm_sur_order$stage == 1), ]
colnames(Rt_Exp_Cli_stageI)[grep("risk_score", colnames(Rt_Exp_Cli_stageI))] <- "risk_score_stageI"
SurMainFun("risk_score_stageI", Rt_Exp_Cli_stageI, DatType =  "ConType", feature = "OS", OutType = "SurF", SurvType = "ALL", len_a = 'high_risk', len_b = "low_risk")

Rt_Exp_Cli_stageII <- rt_mvm_sur_order[which(rt_mvm_sur_order$stage == 2), ]
colnames(Rt_Exp_Cli_stageII)[grep("risk_score", colnames(Rt_Exp_Cli_stageII))] <- "risk_score_stageII"
SurMainFun("risk_score_stageII", Rt_Exp_Cli_stageII, DatType =  "ConType", feature = "OS", OutType = "SurF", SurvType = "ALL", len_a = 'high_risk', len_b = "low_risk")

###########################################mutation analysis###############################################
maf_val = read.maf(maf = "/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/UVM/TCGA_UVM_mutect_somatic.maf")
rt_onco <- maf_val@oncoMatrix
onco_mat <- matrix(unlist(strsplit(colnames(rt_onco), "[.]")), byrow = TRUE, ncol = 7)
colnames(rt_onco) <- paste(onco_mat[, 1], onco_mat[, 2], onco_mat[, 3], sep = "-")
mat_pos <- match(row.names(rt_mvm_sur_order), colnames(rt_onco), nomatch = 0)
rt_onco_order <- rt_onco[, mat_pos]

col = c("Frame_Shift_Del" = color_type[1], "Frame_Shift_Ins" = color_type[8] , "In_Frame_Del" = color_type[2], "In_Frame_Ins" = color_type[3], 
        "Missense_Mutation" = color_type[4], "Nonsense_Mutation" = color_type[5], "Splice_Site" = color_type[6], "Multi_Hit" = color_type[7], 
        "Nonstop_Mutation" = color_type[9], "Translation_Start_Site" = color_type[10])

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  }, 
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Missense_Mutation"], col = NA))
  }, 
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  }, 
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
  }, 
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  }, 
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = col["Multi_Hit"], col = NA))
  }
)

ha1 = HeatmapAnnotation(recurrence = recurrence, gender = gender, stage = stage, 
                        col = list(recurrence= structure(names = unique(recurrence), c('grey', 'black', 'white')),
                                   gender = structure(names = unique(gender), c('grey', 'black')),
                                   stage = structure(names = unique(stage), c('black', 'grey', 'white'))), 
                        os_time =  anno_points(as.numeric(os_time), gp = gpar(col = ifelse( status == 1, "black", "grey")),width = unit(4, "cm"), axis = TRUE), 
                        risk_score = anno_points(as.numeric(risk_score), gp = gpar(col = ifelse( risk_score > median(risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE), 
                        show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE),
                        annotation_height = unit(c(5, 5, 5, 30, 30), "mm"), show_annotation_name = TRUE)
                        

ht <- oncoPrint(rt_onco_order[1:10, ], alter_fun = alter_fun, col = col, column_order = NULL, 
                bottom_annotation = ha1)
pdf(file = paste(cancer, '_snv.pdf', sep = ""), 15, 10)
draw(ht)
dev.off()

###wt and mut boxplot
rt_score_mut <- rbind(rt_onco_order[1:10, ], risk_score)
rt_score_mut[c(grep("Frame", c(rt_score_mut)), grep("Mutation", c(rt_score_mut)), 
               grep("Splice_Site", c(rt_score_mut)), grep("Multi_Hit", c(rt_score_mut)))] = "mut"
rt_score_mut[which(c(rt_score_mut) == "")] = 'wt'

###


#make a boxplot
source("/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/make_boxplot.R")
source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_E.R")
rt_box <- data.frame(t(rt_score_mut), stringsAsFactors = FALSE)
for(gnam in colnames(rt_box)[1:10]){
  rt_box_g <- rt_box[, c(gnam, 'risk_score')]
  colnames(rt_box_g) <- c('group', gnam)
  MakBoxPlot(gnam, 'group', rt_box_g, width = 4, height = 6, theme_E)
}


###make survival base on muttation
dosurvival <- function(snam){
  rt_mvm_sur_order[, snam] <- rt_box[, snam]
  Rt_Exp_Cli <- rt_mvm_sur_order[, c("OS_Status", "OS_Time", snam)]
  diff <- survdiff(Surv(OS_Time, OS_Status) ~Rt_Exp_Cli[, snam], data = Rt_Exp_Cli)
  pValue <- 1-pchisq(diff$chisq, df=1)
  pValueR <- round(pValue, 4)
  fit <- survfit(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli[, snam], data = Rt_Exp_Cli)
  pM_List <- list(Nam = colnames(Rt_Exp_Cli)[3], pValue = pValue, pValue_R = pValueR, fit = fit)
  pValue_R <-  pM_List$pValue_R
  
  pdf(file = paste(snam, "_mut", ".pdf", sep = ""))
  par(mar = c(5, 5, 5, 3))
  plot(pM_List$fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
       main = snam)
  legend(0, 0.17, c('mut', 'wt'), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
  legend("topright", paste("p-value:", pValue_R, sep=" "), cex = 1.5)
  dev.off()
}

snam <- c("GNAQ", "GNA11", "BAP1", "SF3B1", "EIF1AX", "COL14A1", "CYSLTR2")
lapply(snam, dosurvival)






###do pca plot
source('/Users/stead/Documents/SourceTree/R/genomic_personalise_analysis/script/pca_plot.R')
source('/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_F.R')
rt_pca <- rt_mvm_sur_order[, risk_genes_sym]
rt_pca_t <-t(rt_pca)
rt_sam <- data.frame(cbind(row.names(rt_pca), c(rep('low_risk', 40), rep('high_risk', 40)), c(rep('cyan3', 40), rep('red', 40))), stringsAsFactors = FALSE)
colnames(rt_sam) <- c('samples', 'group', 'color')
MakPCA(rt_pca_t, rt_sam = rt_sam, nam = 'risk_score_2', theme_F, type = "pdf", width = NULL, height = NULL)

###get GSEA file 
source('/Users/stead/Documents/SourceTree/R/common_used_codes/tcga_tools/Ensemble2SymbelM.R')
rt_order <- rt_cancer[match(rev(row.names(rt_mvm_sur_order)), row.names(rt_cancer), nomatch = 0), ]
rt_order_sym <- E2SM(t(rt_order))
write.table(rt_order_sym, file = "uvm_riks_score.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

###produce cls file which can be used in GSEA 

cls_1 <- paste(length(rt_order_sym), 2, 1, sep = " ")
cls_2 <- paste(paste("#", 'high_risk', sep = ""), 'low_risk', sep = " ")
cls_3 <- paste(paste(rep('high_risk', 40), collapse = " "), paste(rep('low_risk', 40), collapse = " "), collapse = " ")
cls_file <- rbind(cls_1, cls_2, cls_3)
write.table(cls_file, file = paste('uvm', "group.cls", sep = "_"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

###copy number analysis
rt_uvm_cnv <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/UVM/copy_num/TCGA-UVM.masked_cnv.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rt_uvm_cnv$arm <- NA
rt_uvm_cnv$n.probes <- NA
rt_uvm_cnv_c <- rt_uvm_cnv[, c(1, 2, 6, 3, 4, 7, 5)]
colnames(rt_uvm_cnv_c) <- c("sampleID", "chrom", "arm", "start.pos", "end.pos", 'n.probes', "mean")
sam_cnv <- matrix(unlist(strsplit(rt_uvm_cnv_c$sampleID, '-')), ncol = 4, byrow = TRUE)
rt_uvm_cnv_c$sampleID <- paste(sam_cnv[, 1], sam_cnv[, 2], sam_cnv[, 3], sep = "-")
pos_cnv <- match(rt_uvm_cnv_c$sampleID, row.names(rt_mvm_sur_order), nomatch = 0)
cnv_order <- order(pos_cnv, decreasing = TRUE)
rt_uvm_cnv_o <- rt_uvm_cnv_c[cnv_order, ]

pdf(file = 'cnv_risk.pdf', 8, 10)
plotHeatmap(segments = rt_uvm_cnv_o, upper.lim = 1, lower.lim = -1, chrom = NULL, colors = c('purple3', 'white', 'red'))
dev.off()

###do mut complex heatmap
rt_copy_mut <- cbind(rt_box, rt_mvm_sur_order[, c('OS_Time', "OS_Status")])

BAP1 <- rt_copy_mut$BAP1
SF3B1 <- rt_copy_mut$SF3B1
EIF1AX <- rt_copy_mut$EIF1AX

ha1 = HeatmapAnnotation(BAP1 = BAP1, SF3B1 = SF3B1, EIF1AX = EIF1AX, recurrence = recurrence, 
                        col = list(BAP1= structure(names = unique(BAP1), c('grey', 'black')),
                                   SF3B1 = structure(names = unique(SF3B1), c('black', 'grey')),
                                   EIF1AX = structure(names = unique(EIF1AX), c('grey', 'black')), 
                                    recurrence = structure(names = unique(recurrence), c('grey', 'black', 'white'))), 
                        os_time =  anno_points(as.numeric(os_time), gp = gpar(col = ifelse( status == 1, "black", "grey")),width = unit(4, "cm"), axis = TRUE),
                        risk_score = anno_points(as.numeric(risk_score), gp = gpar(col = ifelse( risk_score > median(risk_score), "red", "blue")),width = unit(4, "cm"), axis = TRUE), 
                        show_legend = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                        annotation_height = unit(c(5, 5, 5, 5, 30, 30), "mm"), show_annotation_name = TRUE)


pdf(file = paste(cancer, "cnv_mut.pdf", sep = ""), 10, 8)
ht_list <- Heatmap(exp, name = "center scaled expression", col = colorRamp2(c(min(exp) , 0, max(exp)), c(" blue", "white", "red")),
                   bottom_annotation = ha1, bottom_annotation_height = unit(10, "cm"), #top_annotation = ha2, top_annotation_height = unit(4, "cm"), 
                   cluster_columns = FALSE, column_dend_reorder = FALSE, 
                   cluster_rows = TRUE, row_dend_reorder = FALSE, 
                   show_row_dend = TRUE, show_column_dend = FALSE,
                   show_row_names = TRUE, show_column_names = FALSE) 
draw(ht_list, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()


