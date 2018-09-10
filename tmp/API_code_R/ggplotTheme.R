ggplotThemeA <- theme_bw()+ 
  theme(plot.title = element_text(size=60, face="bold", hjust = 0.5), 
        strip.text = element_text(size = 60),
        axis.title.x = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 1),
        axis.text.x = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 50, color = "black", vjust = 0.5, hjust = 0.5, angle = 450),
        axis.text.y = element_text(size = 50, color = "black",  vjust = 0.5, hjust = 0.5, angle = 1),
        panel.grid.major = element_blank(), panel.border = element_rect(color='black', size = 1.5),
        panel.grid.minor = element_blank(), axis.ticks = element_line(color = "black", size = 1.5),
        legend.title = element_text(size=40, color="black"),
        legend.key.size = unit(2, "cm"),
        legend.text = element_text(colour = 'black', size = 30))