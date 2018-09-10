
rt <- read.table(file = "/Users/stead/Desktop/GDC/LUAD/merge_result/27_log2exp_pac_heatmap/cancer_T_genes_merge.txt", header = T, row.names = 1, sep = "\t")
gnam <- colnames(rt)[-1]
gnam_SEG <- c("ENSG00000019169.10", "ENSG00000122852.13", "ENSG00000168484.11", "ENSG00000185303.14", 
          "ENSG00000203878.10", "ENSG00000231322.4", "ENSG00000260695.1", "ENSG00000248608.2",
          "ENSG00000168878.15", "ENSG00000112175.7", "ENSG00000166961.13")
rt_SEG <- rt[, c("group", gnam_SEG)]
write.table(rt_SEG, file = "rt_SEG.tsv", col.names = T, row.names = T, sep = "\t")

sort_data <- function(k){
  cnam <-  as.character(unique(rt$group))
  cnam_t <- cnam[c(1, k)]
  data <- rt[c(grep(cnam_t[1], rt$group), grep(cnam_t[2], rt$group)), ]
  data <- data[, c("group", gnam)]
  return(data)
}


count_pvalue <- function(i, k){
  data <- sort_data(k)
  rt_t <- data[, c(1, i)]
  T_Value <- t.test(rt_t[, 2]~group, data = rt_t, paired = FALSE)
  txt <- T_Value$p.value[1]
  return(txt)
  }

                 
count_p_list <- function(k, i){
   i = c(2 :(length(colnames(rt))))
   p_list <- lapply(i, count_pvalue, k)
   data_p <- data.frame(matrix(unlist(p_list), byrow = T, ncol = 1), stringsAsFactors=FALSE)
   rownames(data_p) <- colnames(rt[, gnam])  
   return(data_p)
}

k = c(2 : length(unique(rt$group)))
data_list <- lapply(k, count_p_list, i)
data_p_list <- do.call("cbind", data_list)
colnames(data_p_list) <- as.character(unique(rt$group))[-1]

filter_p <- function(t){
  if(any(data_p_list[t, ] > 0.001) == TRUE ){
    return(NULL)
  } else{
    return(data_p_list[t, ])
  }
}
t <- c(1: length(row.names(data_p_list)))
data_SEG_l <- lapply(t, filter_p)
data_SEG <- do.call("rbind", data_SEG_l)
data_SEG_P <- data_SEG[gnam, ]
write.table(data_SEG_P, file = "data_SEG_P.tsv", row.names = T, col.names = T, sep = "\t")

