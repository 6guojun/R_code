library(wordcloud2)
rt <- read.table(file = "/Users/stead/Desktop/Word.txt", sep = "\t", fill = TRUE, row.names = NULL)

rt$Fre <- rt$Word
rt$Word <- rt$row.names
rt_S <- rt
rt_S <- data.frame(rt_S, stringsAsFactors = FALSE)
row.names(rt_S) <- row.names(demoFreq)[1:length(rt$Word)]
write.table(rt_S, file = "Nordata_Word.txt", row.names = TRUE, col.names = TRUE, sep = "\t")


tiff(filename = "Nordata_Word_Cloud.tiff", width = 512, height = 512)
wordcloud2(rt, size = 1.6)
wordcloud2Output(rt, width = "100%", height = "400px")
dev.off()


letterCloud(demoFreq,"R")

rt <- read.table(file = "/Users/stead/Desktop/Nordata_Word.txt", sep = "\t", fill = TRUE, ,row.names = NULL, stringsAsFactors = FALSE)
rt <- rt[, -1]
row.names(rt) <- row.names(demoFreq)[1:length(rt$Word)]
rt$Fre <- demoFreq$freq
write.table(rt, file = "Nordata_Word.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

letterCloud(rt, word = "R")

library(wordcloud2)
library(webshot)
library(htmlwidgets)
#webshot::install_phantomjs()
rt <- read.table(file = "/Users/stead/Desktop/Nordata_Word.txt", sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
Nordata <- letterCloud(rt, word = "ND", size = 0.3, minRotation = -pi/6, maxRotation = -pi/6,
                       rotateRatio = 1)
#Nordata <- wordcloud2(demoFreq, figPath = "peace.png", size = 1.5, color = "skyblue", backgroundColor="black")
#Nordata <- wordcloud2(rt, size = 1.6)
saveWidget(Nordata,"tmp.html",selfcontained = F)
webshot("tmp.html", "Nordata.pdf", delay =5)


