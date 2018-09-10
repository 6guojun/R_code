library(NMF)
###creat an ExpressionSet with expression data
setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/TCGA")
rt_cli <- read.table(file = "/Users/stead/Desktop/LungCancer_SS/LUAD/cli_data/LUAD_clinicalMatrix", sep = "\t",
                     header = TRUE, row.names = 1 )
rt_exp <- read.table(file = "/Users/stead/Desktop/LungCancer_SS/LUAD/exp_data/LUAD_mRNAmatrix.txt",
                     header = TRUE, row.names = 1, sep = ",")
cli_ele <- c("KRAS", "kras_mutation_found", "kras_mutation_result", "X_OS", "X_OS_IND",  "X_RFS", "X_RFS_IND", "pathologic_M",
             "pathologic_N", "pathologic_T", "pathologic_stage", "number_pack_years_smoked",
             "gender", "age_at_initial_pathologic_diagnosis")
rt_cli_m <- rt_cli[cli_ele]


source("/Users/stead/Documents/SourceTree/R/TCGA_analysis/Match_TCGA_data.R")
rt_cli_p <- GetMatData(rt_exp, rt_cli_m, "01", type = "cli")
rt_exp_p <- GetMatData(rt_exp, rt_cli_m, "01", type = "exp")
Num_I <- c(which(rt_cli_p["pathologic_stage", ] == "Stage I"), which(rt_cli_p["pathologic_stage", ] == "Stage IA"), 
           which(rt_cli_p["pathologic_stage", ] == "Stage IB"))
rt_exp_p_I <- rt_exp_p[, Num_I]
rt_cli_p_I <- rt_cli_p[, Num_I]


express <- rt_exp_p_I[-which(apply(rt_exp_p_I, 1, median) == 0), ]
#colnames(express) <- matrix(unlist(strsplit(colnames(express), "[.]")), ncol = 2, byrow = TRUE)[, 2]
#colnames(express) <- matrix(unlist(strsplit(colnames(express), "[_]")), ncol = 4, byrow = TRUE)[, 4]
genes <- data.frame(rep("A", each = length(express[, 1])), stringsAsFactors = FALSE)
colnames(genes) <- c("Description")
row.names(genes) <- row.names(express)
pheno <- data.frame(colnames(express), stringsAsFactors = FALSE)
pheno$Group[grep("01", colnames(express))] <- "Tumor"
#pheno$Group[grep("N", colnames(express))] <- "Normal"
colnames(pheno) <- c("Sample", "Group")
row.names(pheno) <- pheno$Sample
phenoMetaData <- data.frame(labelDescription=c("Sample name from the file FUSCC_LC_FPKM_ALL.txt", "Group"))
genesMetaData <- data.frame(labelDescription=c("Short description of the gene"))

# create experimental data
expData <- new("MIAME", name = "TCGA_LUAD",
               lab = "TCGA", contact = "shangjunv@163.com",
               title = "LungCancer",
               abstract = "NMF_Package",
               url = "cbio.uct.ac.za",
               other = list(notes = "")
)

es_TCGA_LC <- new('ExpressionSet', exprs = as.matrix(express)
                   , phenoData = new('AnnotatedDataFrame', data = pheno, varMetadata = phenoMetaData)
                   , featureData = new('AnnotatedDataFrame', genes, varMetadata = genesMetaData)
                   , experimentData = expData)

es_TCGA_LC <- es_TCGA_LC[1:50, ]

###NMF analysis
estim.r <- nmf(es_TCGA_LC, 2:6, nrun = 50, seed = 123456)
pdf(file = "tumor_fatorization_rank.pdf")
plot(estim.r)
dev.off()

pdf(file = "tumor_estimation_rank.pdf")
consensusmap(estim.r, annCol = es_FUSCC_LC, labCol = NA, labRow = NA)
dev.off()
V.random <- randomize(es_FUSCC_LC)
estim.r.random <- nmf(V.random, 2:6, nrun = 10, seed = 123456)
plot(estim.r, estim.r.random)