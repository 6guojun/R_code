library(rjson)
library(Biobase)
library(uuid)
library(ggplot2)
JsonFile <- "/Users/stead/Desktop/database_information/TestJson.json"
Jsonlist <- DrawJson(JsonFile)

#complete visualization or statistical analysis
#================================================
MakBox <- function(DataBox, defalut = defalut,  notch = notch, varwidth = varwidth, fill = fill, colour = colour, outlier.colour = outlier, outlier.shape = shape){
  #Make a boxplot
    p = ggplot(data = DataBox, aes(Group,  ExpNam)) +
    xlab(gnam) +  ylab("gene expression")
    if(defalut == TRUE){
      p = p + geom_boxplot()
      print(p)
    } else if(defalut == FALSE){
      p = p + geom_boxplot(notch = notch, varwidth = varwidth, fill = fill, colour = colour, outlier.colour = outlier, outlier.shape = shape)
      print(p)
    }
}


TryCheck <- function(k, KeyBarList){
  #judge whether the Key in our KeyBarList
  #TRUE continue next step, FALSE return a warnming information
  Key <- data.frame(KeyList[k], stringsAsFactors = F)
  if(is.element(colnames(Key), names(KeyBarList)) == TRUE){
    if(is.vector(Key[, 1]) == TRUE){
      if(is.numeric(Key[, 1]) == TRUE){
        if(Key[, 1] < max(KeyBarList[[colnames(Key)]]) | Key[, 1] %in% KeyBarList[[colnames(Key)]] == TRUE){
          NordataWarn <- NA
          OutList <- AnaData@NordataWarn
          return(OutList)
        } else {
          NordataWarn <- paste(colnames(Key), " is not in reference range", sep = "")
          OutList <- AnaData@NordataWarn
          return(OutList)
        }
      } else if (is.character(Key[, 1]) == TRUE){
        if(Key[, 1] %in% KeyBarList[[colnames(Key)]]  == TRUE){
          NordataWarn <- NA
          OutList <- AnaData@NordataWarn
          return(OutList)
        } else {
          NordataWarn <- paste(colnames(Key), " is not in reference range", sep = "")
          OutList <- AnaData@NordataWarn
          return(OutList)
        }
      } else {
        NordataWarn <-  paste(colnames(Key), "  is not a numeric or character", sep = "")
        OutList <- AnaData@NordataWarn
        return(OutList)
      } 
    } else {
      NordataWarn <-  paste(colnames(Key), " is not a vector", sep = "")
      OutList <- AnaData@NordataWarn
      return(OutList)
    }
  } else {
    NordataWarn <- paste(colnames(Key), " is not in KeyBarList", sep = "")
    OutList <- AnaData@NordataWarn
    return(OutList)
  }
}


CheckList <- function(Jsonlist, KeyBarList){
  #check whether the key element of JsonList in KeyBarList
  #JsonList contain the key element
  #KeyBarList is a List defined by author contain the information of all key element
  KeyList <- Jsonlist[KeyEle]
  k <- length(KeyList)
  OutList <- lapply(k, TryCheck, KeyBarList)
    }


MakAnalyze <- function(JsonList, KeyBarList, OutJson){
  #source(RunALlCode.R)
  ResDir <- "/Users/stead/Documents/SourceTree/API_code_R/"
  gnam <- JsonList$gene_ensembl_id
  AnaData <- ConvertData(JsonList)
  DataBox <- GenSortBox(AnaData)
  OutList <- AnaData@GIList
  OutList <- CheckList(JsonList, KeyBarList)
  
  MakBox(DataBox, defalut = JsonList$BoxPat["defalut"],  notch = JsonList$BoxPat["notch"], 
         varwidth = JsonList$BoxPat["varwidth"], fill = JsonList$BoxPat["fill"],
         colour = JsonList$BoxPat["colour"], outlier.colour = JsonList$BoxPat["outlier"], 
         outlier.shape = JsonList$BoxPat["shape"])
  OutList$pLinks <- paste(ResDir, UUIDgenerate(), sep = "")
  return(OutList)
}


