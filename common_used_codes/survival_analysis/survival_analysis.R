### please provide a table containing clinical data and genomic data
### DatType didnot test, so only continuous data were test
### usage: SurMainFun(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"),
###                       feature = c("OS", "DFS"), OutType = c("SurP", "SurF"),
###                       SurvType = SurvType, len_a = len_a, len_b = len_b)
### Jun Shang 
### 2017.07.04
library(survival)

PreSurList <- function(Rt_Exp_Cli, feature = c("OS", "DFS")){
  # Rt_Exp_Cli just contain three columns Time, Status and Gr
  # Ge will be produced by the function of MakSurList
  # Gr were applied to distinguish subgroup
  # if over survival were performed, your names of over survival elements should be OS_Time and OS_Status
  # if DFS were performed, the names of DFS elements shoulbe be DFS_Time and DFS_Status
  if(feature == "OS"){
    Rt_Exp_Cli$OS_Time <- as.numeric(as.character(Rt_Exp_Cli$OS_Time))
    Rt_Exp_Cli$OS_Status <- as.numeric(as.character(Rt_Exp_Cli$OS_Status))
    diff <- survdiff(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli$Gr, data = Rt_Exp_Cli)
    pValue <- 1-pchisq(diff$chisq, df=1)
    pValueR <- round(pValue, 4)
    fit <- survfit(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli$Gr, data = Rt_Exp_Cli)
    pM_List <- list(pValue = pValue, pValue_R = pValueR, fit = fit)
    return(pM_List)
  }
  if(feature == "DFS"){
    Rt_Exp_Cli$DFS_Time <- as.numeric(as.character(Rt_Exp_Cli$DFS_Time))
    Rt_Exp_Cli$DFS_Status <- as.numeric(as.character(Rt_Exp_Cli$DFS_Status))
    diff <- survdiff(Surv(DFS_Time, DFS_Status) ~Rt_Exp_Cli$Gr, data = Rt_Exp_Cli)
    pValue <- 1-pchisq(diff$chisq, df=1)
    pValueR <- round(pValue, 4)
    fit <- survfit(Surv(DFS_Time, DFS_Status) ~ Rt_Exp_Cli$Gr, data = Rt_Exp_Cli)
    pM_List <- list(pValue = pValue, pValue_R = pValueR, fit = fit)
    return(pM_List)
  }
}

MakSurList <- function(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"), feature = c("OS", "DFS")){
  # you can use Gnam to grouping, Gnam can be gene name(continuous variable r logical tracing) or 
  # Mutaion(logical tracing) status or clinical elements(logical tracing)
  
  if(length(which((c("DFS_Time", "DFS_Status", "OS_Time", "OS_Status") %in% colnames(Rt_Exp_Cli))  == TRUE)) >1){
      if(DatType == "ConType" & length(unique(Rt_Exp_Cli[, Gnam])) < 3){
        stop("if yout data type is continuous variable, your unique elements should be over 2")
      } else if(DatType == "ConType" & length(unique(Rt_Exp_Cli[, Gnam] > 3))){
        if(length(which(Rt_Exp_Cli[, Gnam] == median(Rt_Exp_Cli[, Gnam]))) < length(Rt_Exp_Cli[, Gnam])/2){
          Gr <- Rt_Exp_Cli[, Gnam] < median(Rt_Exp_Cli[, Gnam])
          Rt_Exp_Cli$Gr <- Gr
          pM_List <- PreSurList(Rt_Exp_Cli, feature = feature)
          return(pM_List) 
        } else {
          print('over 50% expression is 0')
          return(NULL)
        }
        } else if(DatType == "LogType" & length(unique(Rt_Exp_Cli[, Gnam])) != 2){
        stop("if yout data type is continuous variable, your unique elements should be  2")
          } else if(DatType == "LogType" & length(unique(Rt_Exp_Cli[, Gnam])) == 2){
            Gr <- Rt_Exp_Cli[, Gnam]
            Rt_Exp_Cli$Gr <- Gr
            pM_List <- PreSurList(Rt_Exp_Cli, feature = feature)
            return(pM_List)
            } else {
              return(NULL)
              }
    } else {
      stop("please change your colnames which should contain DFS_Time and DFS_Status or OS_Time and OS_Status")
    }
  }



DrawPvalue <- function(pM_List, Gnam, feature = feature, SurvType = SurvType){
  # SurvType is High mean that the resutls of output is group of high expression with shorter survival time
  # SurvType is Low mean that the resutls of output is group of high expression with longer survival time
  # SurvType is ALL mean that the resutls of output is all
  
  SurvHigh = pM_List$fit$surv[pM_List$fit$strata[1]]
  SurvLow =  pM_List$fit$surv[pM_List$fit$strata[1] + pM_List$fit$strata[2]]
  if(is.null(pM_List)){
    return(NULL)
  }else{
    if(SurvType == "High"){
      if(SurvHigh > SurvLow){
          pValue <- c(Gnam, pM_List$pValue)
          return(pValue)   
        }else{
        return(NULL)
      }
    }else if(SurvType == "Low"){
      if(SurvHigh < SurvLow){
          pValue <- c(Gnam, pM_List$pValue)
          return(pValue)   
        }else{
        return(NULL)
      }
    }else{
        pValue <- c(Gnam, pM_List$pValue)
        return(pValue)   
    }
  }
}


DrawSur <- function(pM_List, Gnam, feature = feature, SurvType = SurvType, len_a = len_a, len_b = len_b){
  # pM_List is a list contain the main element fit and survival analysis pValue
  # SurvType is High mean that the resutls of output is group of high expression with shorter survival time
  # SurvType is Low mean that the resutls of output is group of high expression with longer survival time
  # SurvType is ALL mean that the resutls of output is all
  # Low and High just represent the survival rate of high expression group 
  
  SurvHigh = pM_List$fit$surv[pM_List$fit$strata[1]]
  SurvLow =  pM_List$fit$surv[pM_List$fit$strata[1] + pM_List$fit$strata[2]]
  if(is.null(pM_List)){
    return(NULL)
  }else{
    if(pM_List$pValue < 0.05){
      if(SurvType == "High"){
        if(SurvHigh > SurvLow){
          pValue_R <-  pM_List$pValue_R
          pdf(file = paste(feature, "_", Gnam, ".pdf", sep = ""))
          par(mar = c(5, 5, 5, 3))
          plot(pM_List$fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
               main = Gnam) 
          legend(0, 0.17, c(len_a, len_b), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
          legend("topright", paste("p-value:", pValue_R, sep=" "), cex = 1.5)
          dev.off()
        }else{
          return(NULL)
        }
      }else if(SurvType == "Low"){
        if(SurvHigh < SurvLow){
          pValue_R <-  pM_List$pValue_R
          pdf(file = paste(feature, "_", Gnam, ".pdf", sep = ""))
          par(mar = c(5, 5, 5, 3))
          plot(pM_List$fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
               main = Gnam) 
          legend(0, 0.17, c(len_a, len_b), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
          legend("topright", paste("p-value:", pValue_R, sep=" "), cex = 1.5)
          dev.off()
        }else{
          return(NULL)
        }
      }else{
        pValue_R <-  pM_List$pValue_R
        pdf(file = paste(feature, "_", Gnam, ".pdf", sep = ""))
        par(mar = c(5, 5, 5, 3))
        plot(pM_List$fit, lty = 1:1, lwd = 5, cex.main = 2.5, cex.lab = 2.5, col=c("red","blue"), xlab= ("time (day)"), ylab="surival rate",
             main = Gnam) 
        legend(0, 0.17, c(len_a, len_b), cex = 1.5, lty = 1, lwd = 7, col=c("red","blue"))
        legend("topright", paste("p-value:", pValue_R, sep=" "), cex = 1.5)
        dev.off()
      } 
    }else{
      return(NULL)
    }
  }
}


SurMainFun <- function(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"), 
                       feature = c("OS", "DFS"), OutType = c("SurP", "SurF"), 
                       SurvType = SurvType, len_a = len_a, len_b = len_b){
  # Gnam is the element which may be related with survival 
  # Rt_Exp may the genomic data such expression, methylation or mutation
  # Rt_CLi is the Clinical information data
  # DatTye is the dataformat containing continuous variable r logical tracing
  # SurvType is "High" mean that the resutls of output is group of high expression with shorter survival time
  # SurvType is "Low" mean that the resutls of output is group of high expression with longer survival time
  # SurvType is "ALL" mean that the resutls of output is all
  
  pM_List <- MakSurList(Gnam, Rt_Exp_Cli, DatType = DatType, feature = feature)
  if(OutType == "SurF"){
    DrawSur(pM_List, Gnam, feature = feature, SurvType = SurvType, len_a, len_b)
  } else if(OutType == "SurP"){
    Pvalue <- DrawPvalue(pM_List, Gnam, feature = feature, SurvType = SurvType)
    return(Pvalue)
  } else {
    stop('you should chose SurF or SurP')
  }
}
