library(affy)
library(limma)
library(impute)
library(futile.logger)
library(makecdfenv)
library(gpl15048hursta2a520709cdf)
setwd("/Users/stead/Desktop/LungCancer_SS/affy")
##import phenotype data
phenoData <- read.AnnotatedDataFrame('samples.txt')
pheno <- pData(phenoData)
head (pheno)
#import Annotaion
rt_anno <- read.csv("annotation_affy.csv", head = T, stringsAsFactors = FALSE)
head(rt_anno)
  
  
##RMA normalization
#eset.rma = RMA(Data)
MakAffyExp <- function(pheno, celfile.path, cdfname){
  #rownames(pheno) is the names of all CEL files
  #celfile.path is the path of the GSE*** 
  
  eset_rma <- justRMA(filenames = paste(rownames(pheno), '.CEL', sep = ''), celfile.path = celfile.path, cdfname = cdfname)
  datExpr <- exprs(eset_rma)
  imputed_gene_exp <- impute.knn(datExpr, k = 10, rowmax = 0.5,  colmax = 0.8, maxp = 3000, rng.seed = 362436069)
  datExpr2 <- imputed_gene_exp$data
  return(datExpr2)
}


###datExpr2 <- MakAffyExp(pheno, "./GSE72094", "gpl15048hursta2a520709cdf")
  
PreDesign <- function(Group_A, Group_B, pheno){
  #Group_A is one group such an 'Responder'
  #Group_B is the other group such an 'Nonresponder'
  
  Group <- factor(pheno$group, levels = c(Group_A, Group_B))
  design <- model.matrix(~0 + Group)
  colnames(design) <- c(Group_A, Group_B)
  return(design)
}  


  
datExp <- MakAffyExp(pheno, celfile.path, cdfname)
design <- PreDesign(Group_A, Group_B, pheno)
fit <- lmFit(datExp, design)
contrast_matrix <- makeContrasts(Mut-WT, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, adjust.method = "fdr", p.value = 0.05,
                 lfc = 0 , number = 60000, sort.by = 'logFC')

diff$SYMBEL <- rt_anno$NCBI_ID[match(rownames(diff), rt_anno$Affy_ID)]  

write.table(datExp, file = "Exp.txt", sep = "\t")
write.table(diff, file = "FC_DEG.txt", sep = "\t")
  