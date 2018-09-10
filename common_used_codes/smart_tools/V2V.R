#mat_x contain the second columne which have the element names of vec_y
#usage: V2V(mat_x, vec_y)

V2V <- function(mat_x, vec_y){
  #replace all elements of vec_y with the same elements with other name in mat_x
  #mat_x and vex_y must be vector
  if(!is.null(dim(vec_y)[2])){
    stop('vex_y must be vector')
  } else {
    pos_f = grep(mat_x[1], vec_y)
    
    for(i in 1: dim(mat_x)[1]){
      vec_y <- gsub(mat_x[i, 2], mat_x[i, 1], vec_y)
    }
  }
  return(vec_y)
}
