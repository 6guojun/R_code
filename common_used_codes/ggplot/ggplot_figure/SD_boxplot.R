rt <- read.table(file = "/Users/stead/Desktop/Test/ExpData/CH_FPKM_Log_ab.txt")

rt_A <- rt[, grep("A", colnames(rt))]
rt_B <- rt[, grep("B", colnames(rt))]

count_SD <- function(gnam, rt_data){
    SD_A <- length(colnames(rt_data)) - 1/length(colnames(rt_data))*(sd(rt_data[gnam, ]))
    return(SD_A)
}

gnam_A <- c(row.names(rt_A))
SD_A_List <- lapply(gnam_A, count_SD, rt_A)
SD_A <- c(do.call(rbind, SD_A_List))
rt_SD_A <- data.frame(cbind("A", SD_A), stringsAsFactors = FALSE)
colnames(rt_SD_A) <- c("group", "SD")

gnam_B <- c(row.names(rt_B))
SD_B_List <- lapply(gnam_B, count_SD, rt_B)
SD_B <- c(do.call(rbind, SD_B_List))
rt_SD_B <- data.frame(cbind("B", SD_B), stringsAsFactors = FALSE)
colnames(rt_SD_B) <- c("group", "SD")

rt_SD <- data.frame(rbind(rt_SD_A, rt_SD_B), stringsAsFactors = FALSE)
rt_SD$SD <- as.numeric(rt_SD$SD)

ggplot(rt_SD, aes(fill = group, y = SD, x = group)) + 
  geom_boxplot()
+
  theme_bw() +
  facet_grid(Pnam~Gnam) + 
  theme(strip.text.x = element_text(size=20, face = "bold"),
        strip.text.y = element_text(size=20, face = "bold")) +
  theme(plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 1, face = "bold"),
        axis.text.x = element_text(size = 15, color = "black",  vjust = 0.5, hjust = 0.5, angle = 1, face = "bold"),
        axis.title.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 450, face = "bold"),
        axis.text.y = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 1, face = "bold"),
        legend.title = element_text(colour = 'black', angle = 1, size = 15, hjust = 2, vjust =3, face = "bold"),
        legend.text = element_text(colour = 'black', angle = 1, size = 15, hjust = 2, vjust = 3, face = "bold"),
        panel.grid.major = element_blank(), panel.border = element_rect(),
        axis.ticks = element_line(color = "black", size = 0.5))
dev.off()