library(NMF)
load("object.RData")
# update class definition object <- nmfObject(object)
data(esGolub)

esGolub <- esGolub[1:200, ]
esGolub$Sample <- NULL
estim.r <- nmf(esGolub, 2:6, nrun = 10, seed = 123456)
plot(estim.r)
consensusmap(estim.r, annCol = esGolub, labCol = NA, labRow = NA)
V.random <- randomize(esGolub)
estim.r.random <- nmf(V.random, 2:6, nrun = 10, seed = 123456)
plot(estim.r, estim.r.random)



###creat an ExpressionSet with expression data
outdir <- 'data'
data.dir <- '../data/golub'
express <- read.delim(file.path(data.dir,'ALL_AML_data.txt'), header=FALSE)
genes <- read.delim(file.path(data.dir,'ALL_AML_genes.txt'), header=TRUE, row.names=1)
pheno <- read.delim(file.path(data.dir,'ALL_AML_samples.txt'), header=FALSE, row.names=1)	

# add columns: Name, ALL.AML, Cell
pheno$Sample	<- rownames(pheno)	
pheno$ALL.AML <- as.factor(sub("([a-z]+)_[a-z0-9]+_?.*","\\1",pheno$Sample, ignore.case=TRUE))
pheno$Cell <- as.factor(sub("[a-z]+_[a-z0-9]+_?(.*)","\\1",pheno$Sample, ignore.case=TRUE))	
rownames(express) <- rownames(genes)
colnames(express) <- pheno$Sample
phenoMetaData <- data.frame(labelDescription=c("Sample name from the file ALL_AML_data.txt","ALL/AML status","Cell type"))
genesMetaData <- data.frame(labelDescription=c("Short description of the gene"))

# create experimental data
expData <- new("MIAME", name = "Renaud Gaujoux",
               lab = "UCT NBN", contact = "renaud@cbio.uct.ac.za",
               title = "Golub dataset", 
               abstract = "NMF Package",
               url = "cbio.uct.ac.za",
               other = list(notes = "")
)

esGolub <- new('ExpressionSet', exprs=as.matrix(express)
               , phenoData=new('AnnotatedDataFrame', data=pheno, varMetadata=phenoMetaData)
               , featureData=new('AnnotatedDataFrame', genes, varMetadata=genesMetaData)
               , experimentData = expData)
save(esGolub, file='esGolub.rda')