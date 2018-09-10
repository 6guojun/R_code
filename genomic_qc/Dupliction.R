PathR <- "/Users/stead/Desktop/DL"
setwd(PathR)
SamNam <- "abnormal.txt"
Fnam_Rt <- read.table(file = SamNam)
Fnam <- as.character(Fnam_Rt[, 1])

Draw_Dup <- function(Fnam){
  Snam = unlist(strsplit(Fnam, ".", fixed = TRUE))[1]
  Rnam = unlist(strsplit(Fnam, ".", fixed = TRUE))[2]
  rt <- read.table(file = paste(Fnam, "_Duplication.txt", sep = ""),
                   header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  rt_G <- rt[-c(1, 2, length(rt[, 1])), c(1:3)]
  colnames(rt_G) <- c("level", "duplicated", "toltal")
  rt_G$level <- factor(rt_G$level, levels = unique(rt_G$level))
  rt_G$R <- rep(paste("R", Rnam, sep = ""), length(rt_G[, 1]))
  rt_G$Group <- rep(Snam, length(rt_G[,1]))
  return(rt_G)
}

rt_S_List <- lapply(Fnam, Draw_Dup)
rt_S <- data.frame(do.call(rbind, rt_S_List), stringsAsFactors = FALSE)
rt_S$toltal <- as.numeric(rt_S$toltal)

tiff(filename = paste("Duplication", "_", SamNam, ".tiff", sep = ""), width = 720, height = 512)
ggplot(rt_S, aes(x = level, y = toltal, group = Group, colour = Group)) + 
geom_line() + facet_grid(~R) +  
theme_E
dev.off()

