
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
######correlation 
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.

MakCorPlot <- function(DataCor, Tnam = Tnam, height = 512, width = 512, cex.axis = 2, ...){
  #you should prepare the DataCor which is a dataframe(colnames are names of samples and rownames are names of genes)
  #Tnam is the name of plot, height, width, cex.axi, cex.lab, cex.main and cex.sub are setted to change the figure
    tiff(filename = paste(Tnam,  ".tiff", sep = ""), height = height, width = width)
    p <- pairs(DataCor, upper.panel = panel.cor, diag.panel=panel.hist, cex.axis = cex.axis)#, xlim = c(-15, 15), ylim =c(-15, 15))
    print(p)
    dev.off() 
  }
  



