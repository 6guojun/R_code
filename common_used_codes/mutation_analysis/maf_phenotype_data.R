#usage:

library(dplyr)

ContMafCli <- function(mut_mat, exp_phe_mat, phe_nam){
  mut_mat_mat <- matrix(unlist(strsplit(as.character(mut_mat$Tumor_Sample_Barcode), '-')), ncol = 7, byrow = TRUE)
  mut_mat_cli <- data.frame(cbind(paste(mut_mat_mat[, 1], mut_mat_mat[, 2], mut_mat_mat[, 3], sep = "-"), as.character(mut_mat$Tumor_Sample_Barcode)), stringsAsFactors = FALSE)
  pur_phe_mat <- data.frame(cbind(row.names(exp_phe_mat), exp_phe_mat[, c( "OS_Time", "OS_Status", phe_nam)]), stringsAsFactors = FALSE)
  pur_phe_mat$risk_score[pur_phe_mat[, phe_nam] > median(pur_phe_mat[, phe_nam])] <- 'high_risk'
  pur_phe_mat$risk_score[pur_phe_mat[, phe_nam] <= median(pur_phe_mat[, phe_nam])] <- 'low_risk'
  
  colnames(mut_mat_cli)[1] <- 'Tumor_Sample_Barcode'
  colnames(pur_phe_mat)[1] <- 'Tumor_Sample_Barcode'
  
  mut_mat_cli_phe <- merge(mut_mat_cli, pur_phe_mat, by = 'Tumor_Sample_Barcode')
  mut_mat_cli_phe_d <- mut_mat_cli_phe[, -1]
  colnames(mut_mat_cli_phe_d)[1] <- 'Tumor_Sample_Barcode'
  
  return(mut_mat_cli_phe_d)
}

Addphe2Maf <- function(mut_mat, mat_phe){
  mut_mat_mat <- matrix(unlist(strsplit(as.character(mut_mat$Tumor_Sample_Barcode), '-')), ncol = 7, byrow = TRUE)
  mut_mat_cli <- data.frame(cbind(paste(mut_mat_mat[, 1], mut_mat_mat[, 2], mut_mat_mat[, 3], sep = "-"), as.character(mut_mat$Tumor_Sample_Barcode)), stringsAsFactors = FALSE)
  
  colnames(mut_mat_cli)[1] <- 'Tumor_Sample_Barcode'
  colnames(mat_phe)[1] <- 'Tumor_Sample_Barcode'
  
  mut_mat_phe <- merge(mut_mat_cli, mat_phe, by = 'Tumor_Sample_Barcode')
  mut_mat_phe <- mut_mat_phe[, -1]
  colnames(mut_mat_phe)[1] <- 'Tumor_Sample_Barcode'
  
  return(mut_mat_phe)
}


CountMutPhe <- function(maf_val, maf_phe, clinicalFeatures = clinicalFeatures, gene_list, ggcolor = ggcolor, theme_mut){
  MutPheMat <- data.frame(left_join(left_join(left_join(maf_val@variants.per.sample, maf_val@variant.type.summary, by = 'Tumor_Sample_Barcode'), 
                                              maf_val@variant.classification.summary, by = 'Tumor_Sample_Barcode'),  maf_val@clinical.data, by = 'Tumor_Sample_Barcode'),
                          stringsAsFactors = FALSE)
  
  mut_sta <- c("Variants",  "DEL",  "INS", "SNP", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
               "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site")
  
  for(mnam in mut_sta){
    rt_mut_box <- MutPheMat[, c(mnam, 'risk_score')]
    rt_mut_box <- rt_mut_box[c(which(rt_mut_box$risk_score == "high_risk"), which(rt_mut_box$risk_score == "low_risk")), ]
    
    colnames(rt_mut_box) <- c(mnam, 'group')
    MakBoxPlot(mnam, 'group', rt_mut_box, width = 6, height = 6, theme_mut, ggtype = 'violin', ggcolor = ggcolor, ggylab = 'mutations')
  }
  
  pdf(file = 'mutation_dis_phe.pdf', 10, 10)
  on_ht <- oncoplot(maf = maf_val, clinicalFeatures = clinicalFeatures, sortByAnnotation = TRUE)
  print(on_ht)
  dev.off()
  
  pdf(file = 'sig_genes_mut_dis.pdf', 10, 4)
  sig_genes_ht <- oncostrip(maf = maf_val, genes = gene_list, clinicalFeatures = clinicalFeatures, sortByAnnotation = TRUE)
  print(sig_genes_ht)
  dev.off()
  
  sig_gene_mut <- maf_val@gene.summary$Hugo_Symbol[1:40]
  
  for(gnam in sig_gene_mut){
    pdf(file = paste(gnam, '_mut_sur.pdf', sep = ""))
    sur_mut <- mafSurvival(maf = maf_val, genes = gnam, time = 'OS_Time', Status = 'OS_Status')
    print(sur_mut)
    dev.off()
  }
}

MafToolSum <- function(maf_val, phe_type = c('subtype', 'all'), sub1 = NULL, sub2 = NULL, clinicalFeatures = clinicalFeatures, 
                       gene_list = gene_list, gene_sur = gene_sur, ggcolor = ggcolor, annotationColor = annotationColor, 
                       col_sur = c("maroon", "royalblue")){

  PermafSum <- function(maf_per, P_nam = P_nam){
    pdf(paste('mafSummary_', P_nam, ".pdf", sep = ''))
    plotmafSummary(maf = maf_per, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
  }
  
  Peroncoplot <- function(maf_per, P_nam = P_nam, gene_list = gene_list, clinicalFeatures = clinicalFeatures, annotationColor = annotationColor){
    pdf(paste("oncoplot_", P_nam, ".pdf", sep = ""))
    oncoplot(maf = maf_per, top = 20, fontSize = 12, clinicalFeatures = clinicalFeatures, genes = gene_list, 
             annotationColor = annotationColor)
    dev.off()
  }
  
  Peroncostrip <- function(maf_per, P_nam = P_nam, gene_list = gene_list, clinicalFeatures = clinicalFeatures, annotationColor = annotationColor){
    pdf(paste("oncostrip_", P_nam, ".pdf", sep = ""), 10, 6)
    oncostrip(maf = maf_per, clinicalFeatures = clinicalFeatures, genes = gene_list, annotationColor = annotationColor)
    dev.off()
    }
  
  PerTiTv <- function(maf_per, P_nam = P_nam){
      maf_titv = titv(maf = maf_per, plot = FALSE, useSyn = TRUE)
      #plot titv summary
      pdf(paste("TiTv_", P_nam, ".pdf", sep = ""))
      plotTiTv(res = maf_titv)
      dev.off()
  }
  
  PerVaf <- function(maf_per, P_nam = P_nam, gene_list = gene_list){
    pdf(paste("Vaf_", P_nam, ".pdf", sep = ""))
    plotVaf(maf = maf_per,  genes = gene_list)
    dev.off()
  }
  
  PerInteraction <- function(maf_per, P_nam = P_nam, gene_list = gene_list){
    pdf(paste("Interaction_", P_nam, ".pdf", sep = ""))
    somaticInteractions(maf = maf_per, top = 25, pvalue = c(0.05, 0.1),  genes = gene_list)
    dev.off()
  }
  
  PerSur <- function(gene_sur = gene_sur, maf_per, P_nam = P_nam, col = col_sur){
    pdf(paste(gene_sur, P_nam, "sur.pdf", sep = "_"))
    Sur_F <- mafSurvival(genes = gene_sur, maf = maf_per, time = 'OS_Time', Status = 'OS_Status', isTCGA = FALSE)
    print(Sur_F)
    dev.off()
  }
  
  PerForest <- function(maf_cancer_sub1, maf_cancer_sub2, sub1, sub2){
    pt.vs.rt <- mafCompare(m1 = maf_cancer_sub1, m2 = maf_cancer_sub2, m1Name = sub1, m2Name = sub2, minMut = 5)
    pdf(paste("Forest_", sub1, "_", sub2, ".pdf", sep = ""))
    forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
    dev.off()
  }
  
  PerCoOncoplot <- function(maf_cancer_sub1, maf_cancer_sub2, sub1, sub2, gene_list = gene_list, clinicalFeatures = clinicalFeatures, 
                            annotationColor = annotationColor){
    pdf(paste("CoOncoplot_", sub1, "_", sub2, ".pdf", sep = ""), 10, 6)
    coOncoplot(m1 = maf_cancer_sub1, m2 = maf_cancer_sub2, m1Name = sub1, m2Name = sub2, genes = gene_list, removeNonMutated = TRUE)
    dev.off()
  }
  
  if(phe_type == 'subtype'){
    tsb_sub1 <- maf_val@clinical.data$Tumor_Sample_Barcode[which(maf_val@clinical.data$subtype == sub1)]
    tsb_sub2 <- maf_val@clinical.data$Tumor_Sample_Barcode[which(maf_val@clinical.data$subtype == sub2)]
    
    maf_cancer_sub1 <- subsetMaf(maf_val, tsb = tsb_sub1, restrictTo = "all", mafObj = TRUE)
    maf_cancer_sub2 <- subsetMaf(maf_val, tsb = tsb_sub2, restrictTo = "all", mafObj = TRUE)
    
    for(maf_per in c(maf_val, maf_cancer_sub1, maf_cancer_sub2)){
      P_nam = dim(maf_per@variants.per.sample)[1]
        PermafSum(maf_per, P_nam = P_nam)
        Peroncoplot(maf_per, P_nam = P_nam, gene_list = gene_list, clinicalFeatures = clinicalFeatures, annotationColor = annotationColor)
        Peroncostrip(maf_per, P_nam = P_nam, gene_list = gene_list, clinicalFeatures = clinicalFeatures, annotationColor = annotationColor)
        PerTiTv(maf_per, P_nam = P_nam)
        PerVaf(maf_per, P_nam = P_nam, gene_list = gene_list)
        PerInteraction(maf_per, P_nam = P_nam, gene_list = gene_list)
        lapply(gene_sur, PerSur, maf_per, P_nam = P_nam, col = col_sur)
    }
    
    PerForest(maf_cancer_sub1, maf_cancer_sub2, sub1, sub2)
    PerCoOncoplot(maf_cancer_sub1, maf_cancer_sub2, sub1, sub2, gene_list = gene_list, clinicalFeatures = clinicalFeatures, 
                  annotationColor = annotationColor)
    
  } else if (phe_type == 'all'){
    PermafSum(maf_val, P_nam = "ALL")
    Peroncoplot(maf_val, P_nam = 'ALL', gene_list = gene_list, clinicalFeatures = clinicalFeatures, annotationColor = annotationColor)
    Peroncostrip(maf_val, P_nam = 'ALL', gene_list = gene_list, clinicalFeatures = clinicalFeatures, annotationColor = annotationColor)
    PerTiTv(maf_val, P_nam = 'ALL')
    PerVaf(maf_val, P_nam = 'ALL', gene_list = gene_list)
    PerInteraction(maf_val,  P_nam = 'ALL', gene_list = gene_list)
    lapply(gene_sur, PerSur, maf_per, P_nam = 'ALL', col = col_sur)
  } else {
  print('you must set phe_type')
  }
}

#MafToolSum(maf_sub, phe_type = 'subtype', sub1 = "1", sub2 = "2", clinicalFeatures = NULL, gene_list = NULL, gene_sur = gene_sur, 
#           ggcolor = ggcolor, annotationColor = annotationColor, col_sur = c("maroon", "royalblue"))

