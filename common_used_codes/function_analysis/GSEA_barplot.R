###make a barplot of GSEA results

GSEA_count <-  function(term_type1 = 'gsea_kegg', term_type2 = 'KEGG_', assess_type = c('fdr', 'pval'), output_num = 10, low = 'steelblue1', high = "steelblue4", 
         special_term = NULL, data_dir = data_dir,  width = 10, heith = 5){

  rt_term <- read.table(file = paste(data_dir, term_type1, '/gseapy.gsea.phenotype.report.csv', sep = ''), sep = ",", header = TRUE)
  rt_term$Term <- gsub(term_type2, "", rt_term$Term)
 
   if(is.null(special_term)){
    rt_term = rt_term
  } else {
    rt_term <- rt_term[grep(special_term, rt_term$Term), ]
  }
    
  if(assess_type == 'pval'){
    rt_term_d <- rt_term[which(rt_term$pval < 0.05), c('Term', 'nes', 'pval', 'fdr')]
    
    if(length(rt_term_d[, 1]) > 10){
      rt_term__order <- rt_term_d[order(rt_term_d$pval), ]
      rt_term_num <- rt_term__order[1:output_num, ]
    } else if (length(rt_term_d[, 1]) > 0 & length(rt_term_d[, 1]) < 10) {
      rt_term_num <- rt_term_d
    } else {
      print('no term achieve stander')
    }
    pdf(paste(term_type2, 'pval.pdf'), width, heith)
    pl_bar <- ggplot(rt_term_num, aes(x = reorder(Term, pval), y = nes, fill = pval)) + geom_bar(stat = "identity") +
      scale_fill_gradient(low = low, high = high) + ylab('NES') + xlab(paste(term_type2, 'term', sep = '')) + coord_flip()
    print(pl_bar)
    dev.off()
    
    } else {
      rt_term_d <- rt_term[which(rt_term$fdr < 0.25), c('Term', 'nes', 'pval', 'fdr')]
      
      if(length(rt_term_d[, 1]) > 10){
        rt_term_num <- rt_term_d[1:output_num, ]
      } else if (length(rt_term_d[, 1]) > 0 & length(rt_term_d[, 1]) < 10) {
        rt_term_num <- rt_term_d
      } else {
        print('no term achieve stander')
      }
      
      pdf(paste(term_type2, 'fdr.pdf'), width, heith)
      pl_bar <- ggplot(rt_term_num, aes(x = reorder(Term, fdr), y = nes, fill = fdr)) + geom_bar(stat = "identity") +
        scale_fill_gradient(low = low, high = high) + ylab('NES') + xlab(paste(term_type2, 'term', sep = '')) + coord_flip() 
      print(pl_bar)
      dev.off()
    }
}

data_dir = '/Users/stead/Desktop/GSEA_analysis/new/'

GSEA_count(term_type1 = 'gsea_kegg', term_type2 = 'KEGG_', assess_type = 'fdr', output_num = 19, low = 'steelblue1', high = "steelblue4", 
           special_term = NULL, data_dir = data_dir, width = 10, heith = 5)

GSEA_count(term_type1 = 'gsea_GO', term_type2 = 'GO_', assess_type = 'fdr', output_num = 20, low = 'steelblue1', high = "steelblue4", 
           special_term = NULL, data_dir = data_dir, width = 15, heith = 5)
