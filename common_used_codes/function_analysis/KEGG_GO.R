###you shold prepare a fc and pvalue table which contain gene name, pval and FC 
###you must convert your geneID to ENTREZID before you do kegg and go analysis
###packages:clusterProfiler, org.Hs.eg.db, mygene, annotables, pathview
###ENTREZID: 4312
###ENSEMBLEID: ENSG00000211445
###ENSEMBL: GPX3
###Usemethod
###MakKEGO(rt_KG, ont = ont, data_type = "gene_nam", Gene_ID = "ENSEMBL", OrgDb = "org.Rn.eg.db)
###MakKEGO(rt_KG, ont = ont, data_type = "gene_nam", Gene_ID = "ENSEMBL", OrgDb = "org.Rn.eg.db)
###MakKEGO(rt_KG, ont = ont, data_type = "gene_list", Gene_ID = "SYMBOL", OrgDb = "org.Rn.eg.db)
###species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
###reference: https://bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

library(clusterProfiler)
library(org.Rn.eg.db)
library(mygene)
library(annotables)
library(pathview)

########################################################################
#rt <- read.table(file = "/Users/stead/Desktop/LungCancer_SS/LUAD/DEG/WT_Mut_DEG_Pval_FC.txt", 
#                 header = TRUE, row.names = 1, sep = "\t")
#row.names(rt) <- matrix(unlist(strsplit(row.names(rt), "\\.")), ncol = 2, byrow = TRUE)[, 1]
#######################################################################

##########################KEGG and GO analysis###########################
MakKEGO <- function(rt_KG, ont = ont, data_type = c("gene_nam", "gene_list"),fc_v = fc_v, Gene_ID = c("ENSEMBL", "SYMBOL", "ENTREZID"), OrgDb = OrgDb, organism = organism){
  #if you don't have gene_list which contain FC or log2FC (we need the colname is log2FC), you can set rt_KG as gene names
  #so you will get the results about gene names
  #if you have gene_list wihch contain log2FC, you can get all results
  #ont is one of 'MF', 'BP', 'CC' subontologies 
  
  GetGeneList <- function(rt_KG, data_type = c("gene_nam", "gene_list"), Gene_ID = Gene_ID, fc_v){
    #if you only used gene names as input, your data_type can choice gene_name 
    #the main purpose is convert 'ENSEMBL' or 'SYMBOL' ID to 'ENTREZID'
    #you can choice ENSEMBL or SYMBL based on the type of your gene id 
    #if you produce genelist, youcolnames of rt_KG are gnam and log2FC
    
    if(data_type == "gene_nam"){
      if(Gene_ID == "SYMBOL"){
        gene_df <- bitr(rt_KG, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = OrgDb)
        gene <- gene_df$ENTREZID
        return(gene)
      } else if (Gene_ID == "ENSEMBL"){
        gene_df <- bitr(rt_KG, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = OrgDb)
        gene <- gene_df$ENTREZID
        return(gene)
      } else if (Gene_ID == "ENTREZID") {
        return(rt_KG)
      }
    } else if (data_type == "gene_list"){
      if(Gene_ID == "SYMBOL"){
        gene_df <- bitr(rt_KG$gnam, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = OrgDb)
        gene_list <-  rt_KG$log2FC[match(gene_df$SYMBOL, rt_KG$gnam, nomatch = NA)]
        names(gene_list) <- gene_df$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE)
        return(gene_list)
      } else if (Gene_ID == "ENSEMBL"){
        gene_df <- bitr(rt_KG$gnam, fromType = "ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb = OrgDb)
        gene_list <-  rt_KG$log2FC[match(gene_df$ENSEMBL, rt_KG$gnam, nomatch = NA)]
        names(gene_list) <- gene_df$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE)
        return(gene_list)
      } else if (Gene_ID == "ENTREZID") {
        gene_df <- bitr(rt_KG$gnam, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb = OrgDb)
        gene_list <-  rt_KG$log2FC[match(gene_df$ENTREZID, rt_KG$gnam, nomatch = NA)]
        names(gene_list) <- gene_df$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE)
        return(gene_list)
      } else {
        print ("Gene_ID msut be ENSEMBL, SYMBOL or ENTREZID")
      }
    }
  }
  
###get plot
    if(data_type == "gene_nam"){
      gene = GetGeneList(rt_KG, data_type = "gene_nam", Gene_ID = Gene_ID)
      head(gene)
      eg2np <- bitr_kegg(gene, fromType='ncbi-geneid', toType='kegg', organism= organism)
      
      ggo <- groupGO(gene = gene, OrgDb = OrgDb, ont = ont, level = 3, readable = TRUE)
      write.table(ggo, file = paste(ont,  "_go.txt",  sep = ""), sep = "\t")

      pdf(file = paste(ont, "_ggo",  ".pdf", sep = ""))
      p <- barplot(ggo, drop = TRUE, showCategory = 12)
      print(p)
      dev.off()
      
    } else if (data_type == "gene_list"){
      gene_list =  GetGeneList(rt_KG, data_type = "gene_list", Gene_ID = Gene_ID)
      gene <- names(gene_list)[abs(gene_list) > fc_v]
      head(gene)
      ##
      ggo <- groupGO(gene = gene, OrgDb = OrgDb, ont = ont, level = 3, readable = TRUE)
      head(ggo)
      write.table(ggo, file =  paste(ont, "_go_A.txt", sep = ""), sep = "\t")
      
      pdf(file = paste(ont, "_ggo",  ".pdf", sep = ""))
      p <- barplot(ggo, drop = TRUE, showCategory = 12)
      print(p)
      dev.off()
      ##
      ego <- enrichGO(gene = gene, universe = names(gene_list), OrgDb = OrgDb, ont = ont, pAdjustMethod = "none", 
                      pvalueCutoff  = 0.05, qvalueCutoff  = 0.1, readable = TRUE)
      if(dim(ego)[1] > 0){
        write.table(ego, file =  paste(ont, "_go_B.txt", sep = ""), sep = "\t")
        head(ego)
        pdf(file = paste(ont,  "_ego_barplot.pdf", sep = ""))
        p_A <- barplot(ego, showCategory = 10)
        print(p_A)
        dev.off()
        pdf(file = paste(ont, "_ego_dotplot.pdf", sep = ""))
        p_B <- dotplot(ego)
        print(p_B)
        dev.off()
        pdf(file = paste(ont, "_enrichMap.pdf", sep = ""))
        p_C <- enrichMap(ego)
        print(p_C)
        dev.off()
        pdf(file = paste(ont, "_cnetplot.pdf", sep = ""))
        p_D <- cnetplot(ego, categorySize="pvalue", foldChange = gene_list)
        print(p_D)
        dev.off()
      } else {
       print("enrichGO have no results")
      }
    
      ###KEGG analysis
      kk <- enrichKEGG(gene = gene, organism = organism, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
      write.table(kk, file = "KEGG.txt", sep = "\t")
      pdf(file = paste("kk",  "_barplot.pdf", sep = ""))
      p_A <- barplot(kk, showCategory = 13)
      print(p_A)
      dev.off()
      pdf(file = paste("kk", "_dotplot.pdf", sep = ""))
      p_B <- dotplot(kk, showCategory = 13)
      print(p_B)
      dev.off()
      for(i in kk[, 1]){
        pathview(gene.data  = gene_list,  pathway.id = i, species = organism,
                 limit = list(gene = max(abs(gene_list)), cpd = 1))
      } 
    }
  }
####################################################################

