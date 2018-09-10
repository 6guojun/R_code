##set the Results as the workspace
##you shoud prepare data in the ExpData file
##the files contain FPfile which contain project name, expression data table
##you can choice which project will be shown in figure based on setting PMnam
##Fnam is the expression table name and Pnam is the Project name
##Snam is the samples name(A and B samples of T and N samples)
##before you run this code, you should make the files contained ExpData and script
##SJ 2017/5/17

FPfile <- read.table(file = "./ExpData/FPfile.txt", header = TRUE)
FAnam <- c(as.character(FPfile$Fname))
Pnam <- c(as.character(FPfile$Pname)) ##Pnam can not have replicate
Fnam <- c(FAnam[c(2, 3, 4)])
PMnam <- c(Pnam[c(2, 3, 4)]) ##you can change to choice which Project to shown
Snam <- c("B") ##snam <-c("A", "B", "ab")| c("T", "N", "tn")


source("/Users/stead/Documents/SourceTree/R/genomic_qc/RNA_seq_QC/Correlation/MakCorPlot.R")
source("/Users/stead/Documents/SourceTree/R/genomic_qc/RNA_seq_QC/Correlation/MainMakCor.R")
MainMake(rt_ALL, Fnam, Pnam, PMnam = PMnam,  Snam = Snam)



