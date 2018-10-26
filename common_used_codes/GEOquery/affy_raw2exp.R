setwd('/Users/stead/Downloads/affy_raw_test')
library(affy)
library(limma)
library(impute)
library(diceR)

##import phenotype data
sample_dat = read.table('/Users/stead/Downloads/affy_raw_test/GSE29013_RAW/samples.txt', stringsAsFactors = FALSE, sep = '\t')
sam_nam <- sample_dat$V1

#import Annotaion
anno = read.csv("Annotation.csv", head=T)

##RMA normalization
eset.rma <- justRMA(filenames=paste(sam_nam, '.CEL.gz', sep=''), celfile.path='./GSE29013_RAW', compress = TRUE)
datExpr = exprs(eset.rma)
colnames(datExpr) <- matrix(unlist(strsplit(colnames(datExpr), '_')), ncol = 2, byrow = TRUE)[, 1]
write.table(datExpr, file = '/Users/stead/Desktop/PD-L1_and_TMI_type/GEO_data/NSCLC/GSE29013/expression_data/GSE29013_mat.txt', 
            row.names = TRUE, col.names = TRUE, sep = '\t')
