#Convert the main information of json file to a list
#need library r packages: rjson
----------------------------------------------------
DrawJson <- function(JsonFile){
  RawJson <- fromJSON(file = JsonFile)
  JsonList <- RawJson$gene_exprs
  return(JsonList)
}

