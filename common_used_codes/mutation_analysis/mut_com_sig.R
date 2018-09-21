maf_val = read.maf(maf = "/Users/stead/Desktop/PD-L1_and_TMI_type/UCSC_GDC_data/KIRP/somatic_mut/TCGA.KIRP.mutect.1ab98b62-5863-4440-84f9-3c15d476d523.DR-10.0.somatic.maf")
oncostrip(maf = maf_val, genes = c('TTN', 'KMT2C'))

sam_vars <- data.frame(cbind(maf_val@variants.per.sample, maf_val@variant.classification.summary, maf_val@variant.type.summary), stringsAsFactors = FALSE)
mat_sam_id <- matrix(unlist(strsplit(as.character(sam_vars$Tumor_Sample_Barcode), '-')), ncol = 7, byrow = TRUE)
sam_vars$Tumor_Sample_Barcode <- paste(mat_sam_id[, 1], mat_sam_id[, 2], mat_sam_id[, 3], sep = '-')

mut_score_int <- intersect(sam_vars$Tumor_Sample_Barcode, row.names(rt_mvm_sur_all_min))
mat_mut_score <- data.frame(cbind(sam_vars[match(mut_score_int, sam_vars$Tumor_Sample_Barcode, nomatch = 0), ], 
                                  rt_mvm_sur_all_min$risk_score[match(mut_score_int,  row.names(rt_mvm_sur_all_min), nomatch = 0)]), stringsAsFactors = FALSE) 
colnames(mat_mut_score)[length(mat_mut_score)] <-  'risk_score'
ggplot(mat_mut_score, aes(reorder(Tumor_Sample_Barcode, risk_score), y = SNP)) + geom_bar(stat="identity")

oncoplot(maf_val)


col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

#Color coding for FAB classification; try getAnnotations(x = laml) to see available annotations.
require(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') #path to TCGA LAML MAF file
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') # clinical information containing survival information and histology. This is optional

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)
laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
oncoplot(maf = laml, colors = col, mutsig = laml.mutsig, mutsigQval = 0.01, clinicalFeatures = 'FAB_classification', 
         sortByAnnotation = TRUE, annotationColor = fabcolors)
