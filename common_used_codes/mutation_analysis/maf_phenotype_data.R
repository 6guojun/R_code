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

#mut_cli_phe <- ContMafCli(mut_sam, rt_mvm_sur_all_min)
#maf_val = read.maf(maf = "/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/KIRP/somatic_mut/TCGA.KIRP.mutect.1ab98b62-5863-4440-84f9-3c15d476d523.DR-10.0.somatic.maf", 
#                   clinicalData = mut_cli_phe)

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

#gene_list = c('SETD2', 'PBRM1')
#CountMutPhe(maf_val, mut_cli_phe, clinicalFeatures = 'risk_score', gene_list)



