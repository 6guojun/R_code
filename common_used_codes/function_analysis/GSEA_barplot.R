cancer <- "KIRP"
setwd(paste('/Users/stead/Desktop/subtype_analysis/signature/", cancer, "/function_anlaysis', sep = ""))

kegg_immune <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/immune_related_pathway_and_GO/immune_pathway.txt", sep = "\t",
                          stringsAsFactors = FALSE)
go_immune <- read.table(file = "/Users/stead/Desktop/PD-L1_and_TMI_type/common_used_data/immune_related_pathway_and_GO/go_term.txt", sep = "\t", 
                        stringsAsFactors = FALSE)

rt_kegg <- read.table(file = "GSEA_KEGG.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
rt_kegg_fdr <- rt_kegg[rt_kegg$fdr < 0.25, ]

rt_go <- read.table(file = "GSEA_GO.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
rt_go_fdr <- rt_go[rt_go$fdr < 0.25, ]
rt_go_immune_fdr <- rt_go_fdr[grep("IMMUNE", rt_go_fdr$Term), ]
rt_go_immune_fdr$Term <- gsub("GO_", "", rt_go_immune_fdr$Term)

source("/Users/stead/Documents/SourceTree/R/common_used_codes/ggplot/ggplot_theme/Theme_coord.R")
pdf(file = paste(cancer, "GSEA_GO_immune.pdf", sep = "_"), 20, 10)
ggplot(rt_go_immune_fdr, aes(x = reorder(Term, -fdr) , y = fdr)) + geom_bar(stat="identity", fill="steelblue") + coord_flip() +
   xlab("immune related go term") + ylab("fdr") + theme_coord
dev.off()
