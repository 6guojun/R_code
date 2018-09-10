###make a cox risk modle

MakSum <- function(i, mat_coef, coef_value){
  #i is the length of coef_value
  #the colnames of mat_coef is gene names and the rownames is samples id 
  #coef_value is the coefficient value of each  genes
  
  coef_exp <- mat_coef[, 1]*coef_value[1]
  for(i in 2 : i){
    coef_expi <- mat_coef[, i]*coef_value[i]
    coef_exp <- coef_exp + coef_expi
  }
  return(coef_exp)
}
