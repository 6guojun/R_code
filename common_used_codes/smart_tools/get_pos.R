##usage: GetPos(vec_1, vec_2)
#vec_1 <- c(1, 2, 3)
#vec_2 <- c(rep(1:10, 3), rep(2:10, 3))

GetPos <- function(vec_x, vec_y){
  #attain all positino of elements of vec_x in vec_y
  #vec_x and vex_y must be vector
  if(!is.null(dim(vec_x)[2]) | !is.null(dim(vec_y)[2])){
    stop('vec_x and vex_y must be vector')
  } else {
    pos_f = grep(vec_x[1], vec_y)
    for(i in vec_x){
      pos_i =grep(i, vec_y)
      pos_f = c(pos_f, pos_i)
    }
  }
  return(unique(pos_f))
}


