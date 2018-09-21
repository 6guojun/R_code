GetTumorS <- function(rt_cli){
  mat_sam <- matrix(unlist(strsplit(rt_cli$submitter_id.samples , '-')), ncol = 4, byrow = TRUE)
  T_index <- grep('0..', mat_sam[, 4])
  mat_sam_t <- mat_sam[T_index, ] 
  rt_cli_t <- rt_cli[T_index, ]
  rt_cli_t$submitter_id.samples <- paste(mat_sam_t[, 1], mat_sam_t[, 2], mat_sam_t[, 3], sep = '-')
  return(rt_cli_t)
}
