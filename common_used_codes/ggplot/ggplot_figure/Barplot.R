rt <- read.table(file = "/Users/stead/Desktop/CHQ20170614/NOVO_mapping_rate.txt", 
                  sep = "\t", header  = TRUE, row.names = NULL, stringsAsFactors = FALSE)
colnames(rt) <- c("site", "samples", "total_reads", "mapping_rate", "group")

MakBarPlot <- function(rt_data = rt_data, Vnum = Vnum, Hnum = Hnum, Ftype = c("tiff", "pdf", "jpeg"), 
                       Vnam = NULL, Hnam = NULL, Tnam = NULL, heigth = NULL, width = NULL,
                       themeP = NULL, Fnam = "plot", ClorF = NULL, yaesmin = yaesmin, yaesmax = yaesmax){
  library(ggplot2)
  if(!is.null(Ftype) & !(Ftype %in% c("tiff", "pdf", "jpeg"))){
    stop("Ftype only have tiff, pdf, jepg and null types")
  }
  if(is.null(Ftype)){
    tiff(filename = paste(Fnam, ".tiff", sep = ""), height = heigth, width = width)
    p = ggplot(rt_data, aes(x = get(Vnum), y = get(Hnum), fill = get(ClorF)))  + 
      geom_bar(stat="identity") +  coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) + themeP
    print(p)
    dev.off()  
  }
  if(Ftype == "tiff"){
    tiff(filename = paste(Fnam, ".tiff", sep = ""), height = heigth, width = width)
    p = ggplot(rt_data, aes(x = get(Vnum), y = get(Hnum), fill = get(ClorF)))  + 
      geom_bar(stat="identity") + coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) + themeP
    print(p)
    dev.off()
  }
  if(Ftype == "pdf"){
    pdf(file = paste(Fnam, "pdf", sep = ""), height = heigth, width = width)
    p = ggplot(rt_data, aes(x = get(Vnum), y = get(Hnum), fill = get(ClorF))) +
      geom_bar(stat="identity") + coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) + themeP
    print(p)
    dev.off()
  }
  if(Ftype == "jpeg"){
    pdf(file = paste(Fnam, "jpeg", sep = ""), height = heigth, width = width)
    p = ggplot(rt_data, aes(x = get(Vnum), y = get(Hnum), fill = get(ClorF))) +
      geom_bar(stat="identity") + coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) + themeP
    print(p)
    dev.off()
  }
}
#source("Theme_B.R")
MakBarPlot(rt_data = rt, Vnum = "samples", Hnum = "total_reads", Vnam = "samples", Hnam = "total_reads (%)", 
           Tnam = "mapping_rate", heigth = 480, width = 512, themeP = theme_B, ClorF = "group", yaesmin = 100, 
           yaesmax = 200, Fnam = "Hisat_Novo_total_reads")

