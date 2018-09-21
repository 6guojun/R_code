library(ggplot2)
library(easyGgplot2)

MakBarPlot <- function(mat_ht = mat_ht, Vnum = Vnum, Hnum = Hnum, Vnam = NULL, Hnam = NULL, Tnam = NULL, heigth = NULL, width = NULL,
                       fill_col = "steelblue", themeP = NULL, yaesmin = yaesmin, yaesmax = yaesmax){

    pdf(file = paste(Hnum, ".pdf", sep = ""), height = heigth, width = width)
    p = ggplot(mat_ht, aes(x = get(Vnum), y = get(Hnum))) +
      geom_bar(stat = "identity", fill = fill_col) + coord_cartesian(ylim = c(yaesmin, yaesmax)) +
      labs(x = Vnam, y = Hnam, title = Tnam) +  themeP  
    print(p)
    dev.off()
}







