#draw the data from JsonList which came from DrawJson.R
-------------------------------------------------------
SortBox <- function(x){
  #The main function which was used in generic function
  #Sorting out the data which boxplot need
  DataJson <- data.frame(rbind(cbind("Tumor", unlist(x["tumor_expr_value"])),
                               cbind("Normal", unlist(x["normal_expr_value"]))), stringsAsFactors = F)
  colnames(DataJson) <- c("Group", "ExpNam")
  DataJson[, "ExpNam"] <- as.numeric(DataJson[, "ExpNam"])
  return(DataJson)
}
  
setMethod("GenSortBox", "OCOGArray",
          function(object){
            #draw the main data from class OCOGArray
            DataList <- slot(object, "MData")
            BoxData <- SortBox(DataList)
            return(BoxData)
          }
)