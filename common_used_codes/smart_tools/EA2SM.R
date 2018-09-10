#the tool can be used to convert the expression table with microarray or ensembl id to any type gene id (unique id)
#the duplicated genes will get a mean value or you can choose the gene with max value 
#you can also convert a expression table with duplicated genes to unique with DupAver
###usage: rt_ens_m <- EA2SM(rt_exp_ens, in_type = 'ensembl_gene_id', out_type = 'hgnc_symbol', dup_type = 'max')
###       rt_ens_m <- EA2SM(rt_exp_ens, in_type = 'ensembl_gene_id', out_type = 'hgnc_symbol', dup_type = 'mean')
###       rt_arry_m <- EA2SM(rt_exp_arry, in_type = 'affy_hg_u133_plus_2', out_type = 'hgnc_symbol', dup_type = 'max')
###       rt_arry_m <- EA2SM(rt_exp_arry, in_type = 'affy_hg_u133_plus_2', out_type = 'hgnc_symbol', dup_type = 'mean')

#shangjunv@163.com

library(biomaRt)


###get unique expression table and convert expression value of duplicated genes to mean
DupAver <- function(gnam, rt_sym_dup, dup_type = c('mean', 'max')){
  if(dup_type == 'mean'){
    gvalue <- apply(rt_sym_dup[, -1][grep(gnam, rt_sym_dup$gene_id), ], 2, mean)
    return(gvalue)
  } else if (dup_type == 'max'){
    g_mean <- apply(rt_sym_dup[, -1][grep(gnam, rt_sym_dup$gene_id), ], 1, mean)
    gvalue <- rt_sym_dup[, -1][grep(gnam, rt_sym_dup$gene_id), ][which(g_mean == max(g_mean))[1], ]
    return(gvalue)
  } else {
    stop('you should choose mean or max value when the gene name is not unique')
  }
}
  
###get unqiue gene names table
Dup2UiqMat <- function(rt_sym, dup_type = c('mean', 'max')){
  #get unique genes and duplicated genes
  dup_genes <-  unique(rt_sym$gene_id[duplicated(rt_sym$gene_id)])
  dup_pos_list <- lapply(dup_genes, function(x){grep(x, rt_sym$gene_id)})
  dup_pos <- unlist(dup_pos_list)
  rt_sym_u <-  rt_sym[-dup_pos, -1]
  row.names(rt_sym_u) <- rt_sym$gene_id[-dup_pos]
  rt_sym_dup <- rt_sym[dup_pos, ]
  
  #get unique expression table of rt_sym_dup
  gnam <- as.character(unique(rt_sym_dup$gene_id))
  rt_sym_dup_u_list <- lapply(gnam, DupAver, rt_sym_dup, dup_type = dup_type)
  rt_sym_dup_u <- do.call(rbind, rt_sym_dup_u_list)
  row.names(rt_sym_dup_u) <- gnam
  
  #combine rt_sym_u and rt_sym_dup_u
  rt_m <- data.frame(rbind(rt_sym_u, rt_sym_dup_u), stringsAsFactors = FALSE)
  return(rt_m)
}

EA2SM <- function(rt_exp, in_type = in_type, out_type = out_type, dup_type = c('mean', 'max')){
  print('in_type/out_type contain: ensembl_gene_id, affy_hg_u133_plus_2, entrezgene, hgnc_symbol')
  
  if(length(matrix(unlist(strsplit(row.names(rt_exp), "[.]"))))/length(row.names(rt_exp)) == 2){
      row.names(rt_exp) <- matrix(unlist(strsplit(row.names(rt_exp), "[.]")), ncol = 2, byrow = TRUE)[, 1]
      print('esembl id with point')
    } else if (length(matrix(unlist(strsplit(row.names(rt_exp), "[.]"))))/length(row.names(rt_exp)) < 2) {
      row.names(rt_exp) <- row.names(rt_exp)
      print('esembl id without point or other gene id')
    } else {
      stop('check your gene id')
    }
  
    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    eg_arry <- getBM(attributes = c(in_type, out_type), values = row.names(rt_exp), mart = ensembl)
    eg_arry_d <- eg_arry[-c(which(eg_arry[, in_type] == ""), which(eg_arry[, out_type] == "")),  ]
    ens_int <- intersect(eg_arry_d[, in_type], row.names(rt_exp))
    rt_sym <- data.frame(cbind(eg_arry_d[, out_type][match(ens_int, eg_arry_d[, in_type], nomatch = 0)], 
                               rt_exp[match(ens_int, row.names(rt_exp), nomatch = 0), ]), stringsAsFactors = FALSE)
    colnames(rt_sym)[1] <- 'gene_id'
    rt_m <- Dup2UiqMat(rt_sym, dup_type = dup_type)
    return(rt_m)
  }

