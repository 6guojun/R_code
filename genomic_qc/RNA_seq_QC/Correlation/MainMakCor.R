##the main funtiong of the code is make a correlation plot
##you should prapare FPfile which contain the table name and project name
##you can choice which project and which samples to be sent out by setting parameters: Fnam, Pnam, PMnam,  Snam1, Snam2
##Fnam represent the file name, Pnam represent Project name
##PMnam is the main project names to show in correlation figure
##Snam reprsent kinds of samples(the kinds of samples should not over two)
##SJ 2017/05/16
MergeData <- function(Fnam = Fnam, Pnam = Pnam){
  ##when you have many table to merge, you can merge all table by first column
  ##Fnam is the filename, Pname is the Project name
  if(is.null(Fnam) | is.null(Pnam)){
    stop("Fnam or Pnam can not be NULL")
  }
  if(length(Fnam) == 1 ){
    k <- length(Fnam)
      data <- read.table(file = paste("./ExpData/", Fnam, ".txt", sep = ""), header = TRUE,
                           row.names = NULL, stringsAsFactors = FALSE, sep = "\t")
      return(rt_data) 
    }  
  if(length(Fnam) > 1){
    k <- length(Fnam)
    rt_data <- read.table(file = paste("./ExpData/", Fnam[1], ".txt", sep = ""), header = TRUE,
                          row.names = NULL, stringsAsFactors = FALSE, sep = "\t")
    for(i in 2: k){
      x <- Fnam[i]
      y <- Pnam[i]
      data <- read.table(file = paste("./ExpData/", x, ".txt", sep = ""), header = TRUE,
                         row.names = NULL, stringsAsFactors = FALSE, sep = "\t")
      rt_data <- merge(rt_data, data, by =1)
    }
    rt_data <- data.frame(rt_data[, -1], stringsAsFactors = FALSE)
    return(rt_data) 
  }
}

rt_all <- MergeData(Fnam, Pnam)

DrawData <- function(rt_ALL, PMnam = PMnam, Snam = Snam){
  #rt_ALL is the dataframe
  #PMnam is Project name which your choice
  #Snam is the sample name
  #you can choice which project data you want to get based on PMnam
  #you can choice which samples you want to get based on Snam
  if(is.null(PMnam) | is.null(Snam)){
    #you can not set PMnam as NUll
    stop("PCnam, Snam can not be NULL")
  }
  if(length(PMnam) == 1){
    rt_P <- rt_ALL[, grep(PMnam, colnames(rt_ALL))]
    if(length(Snam) == 1){
      rt_P_S <- rt_P[, grep(Snam, colnames(rt_P))]
      return(rt_P_S)
    }
    if(length(Snam) > 1){
      k <- length(Snam)
      rt_P_S <- rt_P[, grep(PMnam[1], colnames(rt_P))]
      for(i in 2: k){
        rt_P_K <- rt_P[, grep(PMnam[i], colnames(rt_P))]
        rt_P_S <- cbind(rt_P_S, rt_P_K)
      }
      rt_P_S <- data.frame(rt_P_S, stringsAsFactors = FALSE)
      return(rt_P_S)
    }
  }
  if(length(PMnam) > 1){
    rt_P <- rt_ALL[, grep(PMnam[1], colnames(rt_ALL))]
    t <- length(PMnam)
    for(r in 2 :t){
      rt_P_T <- rt_ALL[, grep(PMnam[r], colnames(rt_ALL))]
      rt_P <- cbind(rt_P, rt_P_T)
    }
    rt_P <- data.frame(rt_P, stringsAsFactors = FALSE)
    if(length(Snam) == 1){
      rt_P_S <- rt_P[, grep(Snam, colnames(rt_P))]
      return(rt_P_S)
    }
    if(length(Snam) > 1){
      v <- length(Snam)
      rt_P_S <- rt_P[, grep(Snam[1], colnames(rt_P))]
      for(y in 2 :V){
        rt_P_V <- rt_P[, grep(Snam[y], colnames(rt_P))]
        rt_P_S <- cbind(rt_P_S, rt_P_V)
      }
      rt_P_S <- data.frame(rt_P_S, stringsAsFactors = FALSE)
      return(rt_P_S)
    }
  }
  }

rt_Data_C <- DrawData(rt_all, PMnam = PMnam, Snam = Snam)



MainMake <- function(rt_ALL, Fnam = Fnam, Pnam = Pnam, PMnam = PMnam,  Snam = Snam){
  rt_ALL <- MergeData(Fnam, Pnam)
  rt_Data_M <- DrawData(rt_ALL, PMnam, Snam)
  k = length(PMnam)
  Mnam <- PMnam[1]
  for(i in 2: k){
    Mnam <- paste(Mnam, "_", PMnam[i], sep = "")
  }
  Tnam <- paste(Mnam, Snam, sep = "_")
  MakCorPlot(rt_Data_M, Tnam = Tnam)
}
