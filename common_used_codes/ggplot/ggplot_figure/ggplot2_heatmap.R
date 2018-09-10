library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")



mcp_score_mat_gg <- melt(t(mcp_score_mat))
pdf(file = "ggplot2_test.pdf", 20, 5)
heatmap_plot <- ggplot(data = mcp_score_mat_gg, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  theme(axis.text.y = element_text(size = 6)) + theme(plot.title = element_text(size = 3, hjust = 0.5),
                                                      axis.title.x = element_text(size = 3, color = "black", angle = 1),
                                                      axis.text.x = element_text(size = 3, color = "black", angle = 90),
                                                      axis.title.y = element_text(size = 20, color = "black", angle = 1),
                                                      axis.text.y = element_text(size = 20, color = "black", angle = 1))
print(heatmap_plot)
dev.off()


ggplot(rt_mcp_socre_risk, aes(x = group, y = `Monocytic lineage`)) + geom_boxplot()



