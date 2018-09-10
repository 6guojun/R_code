library(NMF)
setwd("/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp")
###creat an ExpressionSet with expression data
DEG_FC <- read.table(file = "/Users/stead/Desktop/FUSCC_LC_QC/FUSCC_LC_RNAseq_ALL/exp/ALL_DEG/Tumor_Normal_DEG_Pval_FC.txt", sep = "\t")

###top100 gene list
DEG_FC_Order <- DEG_FC[order(DEG_FC$log2FC, decreasing = FALSE), ]
Top_genes <- c(row.names(DEG_FC_Order)[1 : 50], row.names(DEG_FC_Order)[(length(DEG_FC_Order$Pval)-50) :  length(DEG_FC_Order$Pval)])

###Top100 express
###get tumor group
express <- read.table(file = "FUSCC_LC_FPKM_ALL.txt", header = TRUE, row.names = 1, sep = "\t")
express <- express[, -c(405:413)]
express <- express[-which(apply(express, 1, median) == 0), ]
express <- express[Top_genes, ]
colnames(express) <- matrix(unlist(strsplit(colnames(express), "[.]")), ncol = 2, byrow = TRUE)[, 2]
colnames(express) <- matrix(unlist(strsplit(colnames(express), "[_]")), ncol = 4, byrow = TRUE)[, 4]
express <- express[, grep("T", colnames(express))]

###make expresseionset object
genes <- data.frame(rep("A", each = length(express[, 1])), stringsAsFactors = FALSE)
colnames(genes) <- c("Description")
row.names(genes) <- row.names(express)
pheno <- data.frame(colnames(express), stringsAsFactors = FALSE)
pheno$Group[grep("T", colnames(express))] <- "Tumor"
pheno$Group[grep("N", colnames(express))] <- "Normal"
colnames(pheno) <- c("Sample", "Group")
row.names(pheno) <- pheno$Sample
phenoMetaData <- data.frame(labelDescription=c("Sample name from the file FUSCC_LC_FPKM_ALL.txt", "Group"))
genesMetaData <- data.frame(labelDescription=c("Short description of the gene"))

# create experimental data
expData <- new("MIAME", name = "FUSCC_PGx",
               lab = "PGx", contact = "shangjunv@163.com",
               title = "LungCancer", 
               abstract = "NMF_Package",
               url = "cbio.uct.ac.za",
               other = list(notes = "")
)

es_FUSCC_LC <- new('ExpressionSet', exprs = as.matrix(express)
               , phenoData = new('AnnotatedDataFrame', data = pheno, varMetadata = phenoMetaData)
               , featureData = new('AnnotatedDataFrame', genes, varMetadata = genesMetaData)
               , experimentData = expData)


###NMF analysis
estim.r <- nmf(es_FUSCC_LC, 2:7, nrun = 50, seed = 123456)
pdf(file = "estim_r.pdf")
plot(estim.r)
consensusmap(estim.r, annCol = es_FUSCC_LC, labCol = 1, labRow = 1)
dev.off()

#estim.r <- nmf(es_FUSCC_LC, 2:3, nrun = 50, seed = 123456)
#plot(estim.r)
#pdf(file = "estim_r.pdf")
#consensusmap(estim.r, annCol = es_FUSCC_LC, labCol = 1, labRow = 1)
#dev.off()

V.random <- randomize(es_FUSCC_LC)

estim.r.random <- nmf(V.random, 2:6, nrun = 10, seed = 123456)
plot(estim.r, estim.r.random)


###visualization methods
res <- nmf(es_FUSCC_LC, 4) #.options = "t"
plot(res)
plot(res.nulti.method)

layout(cbind(1, 2))
#basismap
pdf(file = "basismap.pdf")
basismap(res, Rowv = TRUE)
dev.off()

#coefmap
pdf(file = "coefmap.pdf")
coefmap(res, Colv = FALSE, Rowv = FALSE, tracks = "basis")
dev.off()
