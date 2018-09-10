###Somatic mutation analysis
CountGeneMutRate <- function(mut_mat, mtype = c('snv', 'cnv')){
  #mut_mat:colnames is gene name and row.names is samples id 
  #group contian two column contain samples id and group
  #contain risk
  mut_mat_low_risk <- mut_mat[which(mut_mat$risk_score <= median(mut_mat$risk_score)), ]
  mut_mat_high_risk <- mut_mat[which(mut_mat$risk_score > median(mut_mat$risk_score)), ]
  
  
  GetAllWTPos <- function(rp, mut_mat_g, mtype){
    if(mtype == 'cnv'){
      res_i_1 = 0
      for(i in rp){
        if(any(mut_mat_g[, i] == 'cnv')){
          res_i_1 = res_i_1 + 1
        }else {
          next
        }
      }
      return(res_i_1)
    }
    if(mtype == 'snv'){
      res_i_1 = 0
      for(i in rp){
        if(any(mut_mat_g[, i] == 'mut')){
          res_i_1 = res_i_1 + 1
        }else {
          next
        }
      }
      return(res_i_1)
    }
  }
  
  #get low risk mutation genes
  rp_low = 1:length(mut_mat[1, ])
  mut_pos_low <- GetAllWTPos(rp_low, mut_mat_low_risk, mtype)
  
  #get high risk mutation genes
  rp_high = 1:length(mut_mat[1, ])
  mut_pos_high <- GetAllWTPos(rp_high, mut_mat_high_risk, mtype)
  
  ###get rate
  mut_low_rate <- (mut_pos_low - 1)/length(mut_mat_low_risk[, 1])
  mut_high_rate <- (mut_pos_high - 1)/length(mut_mat_high_risk[, 1])
  
  ###get mut mat
  mat_mut_rate <- data.frame(c(mut_low_rate, mut_high_rate), c('low_risk', 'high_risk'))
  colnames(mat_mut_rate) <- c('gene_mutation_rate', 'group')
  
  pdf(file = paste(cancer, "gene_risk_", mtype, "rate.pdf", sep = "_"), 5, 5)
  mut_bar <- ggplot(mat_mut_rate, aes(x = group, y = gene_mutation_rate)) + geom_bar(stat="identity", fill = c('red', 'blue')) + 
    theme_B
  print(mut_bar)
  dev.off()
}




