theme_cor_heatmap <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(size = 30, colour = 'black', vjust = 0.5, hjust = 0.5, angle = 90), 
  axis.title.y = element_blank(),
  axis.text.y = element_text(size = 30, colour = 'black', vjust = 0.5, hjust = 0.5, angle = 1), 
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.text = element_text(colour = 'black', angle = 1, size = 20, hjust = 2, vjust = 3), 
  legend.position = c(1.19, 1.1),
  legend.direction = "horizontal")
