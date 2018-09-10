#Convert the JsonList to an object contain the class information
#================================================================
ConvertData <- function(JsonList){
  #draw the main data from json file and convert the main data to a class 
  #GIList contain the specific information of MData
  #Base on ClassName of JsonList to choose class
  #AnaData is an object which were return
  GIField = Jsonlist[c("species", "project_name", "tissue_type", "tissue_t_full_name", "gene_ensembl_id")]
  if(JsonList$ClassName == "OCOGArray"){
    AnaData <- new("OCOGArray", MData = JsonList, GIList = GIField)
    return(AnaData) 
  } else if (JsonList$ClassName == "OCMGArray"){
    AnaData <- new("OCMGArray", MData = JsonList, GIList = GIField)
    return(AnaData) 
  } else if (JsonList$ClassName == "MCMGArray"){
    AnaData <- new("MCMGArray", MData = JsonList, GIList = GIField)
    return(AnaData) 
  }
}

AnaData <- ConvertData(JsonList)
