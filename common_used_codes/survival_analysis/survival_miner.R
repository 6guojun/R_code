###do survival analysis 
###usage:source("/Users/stead/Documents/SourceTree/R/common_used_codes/survival_analysis/ggsurv_theme/ggsurv_theme_A.R")
####      DrawSurminer("sex", Rt_Exp_Cli, DatType = "LogType", theme_surv_a, len_a = 'len_a', len_b = 'len_b', palette = paletterisk.table = risk.table)
library(survminer)
library(survival)

#########################test###################################
#setwd("/Users/stead/Desktop/subtype_analysis/survival_test")
#require("survival")
#Rt_Exp_Cli <- lung
#colnames(Rt_Exp_Cli)[c(2, 3)]<- c("OS_Time", "OS_Status")
################################################################


DrawSurminer <- function(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"), theme_surv, len_a = len_a, len_b = len_b, s_cutoff = NULL, 
                         palette = palette, risk.table = risk.table, fun = NULL){
  # pM_List is a list contain the main element fit and survival analysis pValue
  # SurvType is High mean that the resutls of output is group of high expression with shorter survival time
  # SurvType is Low mean that the resutls of output is group of high expression with longer survival time
  # SurvType is ALL mean that the resutls of output is all
  # Low and High just represent the survival rate of high expression group 
  print("your survival data title must be OS_Status and OS_Time")
  print("len_a < median")
  Gnam <<- Gnam
  Rt_Exp_Cli <<- Rt_Exp_Cli
  
  PreSurList <- function(Gnam, Rt_Exp_Cli){
    # Rt_Exp_Cli just contain three columns Time, Status and Gr
    # Ge will be produced by the function of MakSurList
    # Gr were applied to distinguish subgroup
    # if over survival were performed, your names of over survival elements should be OS_Time and OS_Status
    # if DFS were performed, the names of DFS elements shoulbe be DFS_Time and DFS_Status
    
    Rt_Exp_Cli$OS_Time <- as.numeric(as.character(Rt_Exp_Cli$OS_Time))
    Rt_Exp_Cli$OS_Status <- as.numeric(as.character(Rt_Exp_Cli$OS_Status))
    diff <- survdiff(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli[, Gnam], data = Rt_Exp_Cli)
    pValue <- 1-pchisq(diff$chisq, df=1)
    fit <- survfit(Surv(OS_Time, OS_Status) ~ Rt_Exp_Cli[, Gnam], data = Rt_Exp_Cli)
    pM_List <- list(pValue = pValue, fit = fit)
    return(pM_List)
  }
  
  MakSurList <- function(Gnam, Rt_Exp_Cli, DatType = c("ConType", "LogType"), len_a = len_a, len_b = len_b, s_cutoff = NULL){
    # you can use Gnam to grouping, Gnam can be gene name(continuous variable r logical tracing) or 
    # Mutaion(logical tracing) status or clinical elements(logical tracing)
    
    if(length(which((c("DFS_Time", "DFS_Status", "OS_Time", "OS_Status") %in% colnames(Rt_Exp_Cli))  == TRUE)) >1){
        if(DatType == "ConType" & length(unique(Rt_Exp_Cli[, Gnam])) < 3){
          stop("if yout data type is continuous variable, your unique elements should be over 2")
        }
        if(DatType == "ConType" & length(unique(Rt_Exp_Cli[, Gnam])) > 3){
          if(is.null(s_cutoff)){
            Rt_Exp_Cli[, Gnam][Rt_Exp_Cli[, Gnam] > median(Rt_Exp_Cli[, Gnam])] <- len_b
            Rt_Exp_Cli[, Gnam][which(Rt_Exp_Cli[, Gnam] != len_b)] <- len_a
            pM_List <- PreSurList(Gnam, Rt_Exp_Cli)
            return(pM_List)  
          } else if (!is.null(s_cutoff)){
            Rt_Exp_Cli[, Gnam][Rt_Exp_Cli[, Gnam] > s_cutoff]  <- len_b
            Rt_Exp_Cli[, Gnam][which(Rt_Exp_Cli[, Gnam] != len_b)] <- len_a
            pM_List <- PreSurList(Gnam, Rt_Exp_Cli)
            return(pM_List)  
          } 
        }
      if(DatType == "LogType" & length(unique(Rt_Exp_Cli[, Gnam])) < 2){
          stop("if yout data type is continuous variable, your unique elements should be over  2")
        } else if(DatType == "LogType" & length(unique(Rt_Exp_Cli[, Gnam])) >= 2){
          pM_List <- PreSurList(Gnam, Rt_Exp_Cli)
          return(pM_List)
        }
    }else{
      stop("please change your colnames which should contain DFS_Time and DFS_Status or OS_Time and OS_Status")
    }
  }
  
  pM_List <- MakSurList(Gnam, Rt_Exp_Cli, DatType, len_a = len_a, len_b = len_b, s_cutoff = s_cutoff)
  fit <- pM_List$fit
  print(fit)
  pValueR <- pM_List$pValue
  print(names(fit$strata))
  names(fit$strata) <- gsub("Rt_Exp_Cli., Gnam]=", "", names(fit$strata))
  
  if(is.null(fun)){
  pdf(file = paste(Gnam, ".pdf", sep = ""), 11, 11, onefile = FALSE)
   dt <- ggsurvplot(fit, title = Gnam, pval = paste('P = ', format(pValueR, digits = 4, scientific = TRUE), sep = ""), 
                    conf.int = FALSE, pval.size = 10, pval.coord = c(1, 0), xlab = "Time in days", 
                    size = 2, #change line size
                    risk.table = risk.table, # Add risk table
                    palette = palette, # custom color palette
                    risk.table.col = "strata", # Change risk table color by groups
                    ggtheme = theme_surv, fontsize = 8, risk.table.height= 0.3, fun = fun)
   print(dt)
   dev.off()
   
 } else {
 pdf(file = paste(Gnam, ".pdf", sep = ""), 11, 11, onefile = FALSE)
  dt <- ggsurvplot(fit,  title = Gnam, pval = paste('P = ', format(pValueR, digits = 4, scientific = TRUE), sep = ""), 
                   conf.int = FALSE, pval.size = 10, pval.coord = c(1, 1), xlab = "Time in days", 
                    size = 2, #change line size
                    risk.table = risk.table, # Add risk table
                    palette = palette, # custom color palette
                    risk.table.col = "strata", # Change risk table color by groups
                    ggtheme = theme_surv, fontsize = 8, risk.table.height = 0.3, fun = fun)
   print(dt)
   dev.off()
 }
  
  rm(Gnam, Rt_Exp_Cli, pos = ".GlobalEnv")
}

###stratified survival analysis such such as TNM stage

StratifySur <- function(rt_cli_mat, sur_nam = sur_nam, pheno_ele = pheno_ele, strat_ele = strat_ele, title_nam = title_nam, theme_surv, 
                        DatType = c('ConType', 'LogType'), len_a = len_a, len_b = len_b, palette = c('red', 'blue'), risk.table = TRUE, fun = NULL){
  
  print("pheno_ele like stage/TMN, sur_nam like gene or risk score, title_nam like risk_score_M0")
  
  rt_cli_mat_strat <- rt_cli_mat[which(rt_cli_mat[, pheno_ele] == strat_ele), ]
  colnames(rt_cli_mat_strat)[which(colnames(rt_cli_mat_strat) == sur_nam)] <- title_nam
  
  str_sur <- DrawSurminer(title_nam, rt_cli_mat_strat, DatType = DatType, theme_surv, len_a = len_a, len_b = len_b, 
                          s_cutoff = NULL,  palette = palette, risk.table = risk.table, fun = fun)
}
