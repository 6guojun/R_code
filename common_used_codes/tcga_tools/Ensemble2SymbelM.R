#the tool can be used to convert the expression table with ensemble id to gene symbol 
#the duplicated genes will get a mean value
#the tool is designed for TCGA expression table of which colnames is samples id and row.names is ensemble id 
###usage: rt_u <- E2SM(rt_E,  ens_type = c('WP', "NP"))
#you can also convert a expression table with duplicated genes to unique with DupAver
###usage: gnam <- unique(rt_sym_dup$SYMBOL)
#         rt_sym_dup_u_list <- lapply(gnam, DupAver, rt_sym_dup)
#         rt_sym_dup_u <- do.call(rbind, rt_sym_dup_u_list)
#         row.names(rt_sym_dup_u) <- gnam
#shangjunv@163.com


library(clusterProfiler)

###get unique expression table and convert expression value of duplicated genes to mean
DupAver <- function(gnam, rt_sym_dup){
  gvalue <- apply(rt_sym_dup[, -1][grep(gnam, rt_sym_dup$SYMBOL), ], 2, mean)
  return(gvalue)
}


E2SM <- function(rt_E, ens_type = c('WP', "NP")){
  if(ens_type == "WP"){
    ensemble_gene <- matrix(unlist(strsplit(row.names(rt_E), "[.]")), ncol = 2, byrow = TRUE)[, 1]
  } else if (ens_type == "NP") {
    ensemble_gene <- row.names(rt_E)
  }
    
  eg = bitr(ensemble_gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = FALSE)
  ens_m <- intersect(ensemble_gene, eg$ENSEMBL)
  ens_pos <- match(ens_m, ensemble_gene, nomatch = 0)
  sym_pos <- match(ens_m, eg$ENSEMBL, nomatch = 0)
  sym_gene <- eg$SYMBOL[sym_pos]
  
  #remove nomatch ensemble genes
  rt_sym <- data.frame(cbind(sym_gene[!is.na(sym_gene)], rt_E[ens_pos, ][!is.na(sym_gene), ]), stringsAsFactors = FALSE)
  colnames(rt_sym)[1] <- "SYMBOL"
  
  #get unique genes and duplicated genes
  dup_genes <-  unique(rt_sym$SYMBOL[duplicated(rt_sym$SYMBOL)])
  dup_pos_list <- lapply(dup_genes, function(x){grep(x, rt_sym$SYMBOL)})
  dup_pos <- unlist(dup_pos_list)
  rt_sym_u <-  rt_sym[-dup_pos, -1]
  row.names(rt_sym_u) <- rt_sym$SYMBOL[-dup_pos]
  rt_sym_dup <- rt_sym[dup_pos, ]
  
  #get unique expression table of rt_sym_dup
  gnam <- unique(rt_sym_dup$SYMBOL)
  rt_sym_dup_u_list <- lapply(gnam, DupAver, rt_sym_dup)
  rt_sym_dup_u <- do.call(rbind, rt_sym_dup_u_list)
  row.names(rt_sym_dup_u) <- gnam

  #combine rt_sym_u and rt_sym_dup_u
  rt_m <- data.frame(rbind(rt_sym_u, rt_sym_dup_u), stringsAsFactors = FALSE)
}

