#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step06-05-05_tabulate-MRPRESSO-results.R
# Modified from : MR_step06-03-05_tabulate-two-sample-MR-analysis-results.R
# Date created  : 20190420
# Purpose       : Combine MR-PRESSO analysis results as single files
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency:     
# Function external:  
# Type  Files
#----------------------------------------------------------------------------------------
# Input Sys.glob(paste0(loc.twoSampleMR.output,"MR-analysis_exposure-*_outcome-*.tsv")) # 180 files
# Outpu paste0(loc.MRPRESSO.tabulated,"MR-PRESSO-global-test-results",".tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190822  Added 1 more trait pair. Exported paste0(loc.MRPRESSO.tabulated,"MR-PRESSO-global-test-results",".tsv")
# 20190821 Converted causal estimates to represent a doubling in odds of binary exposure, exported paste0(loc.MRPRESSO.tabulated,"MR-PRESSO-global-test_continous-outcomes",".tsv")
# 20190818  Exported the 1 file above  
# 20190420  Exported the 1 file above
#----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Folder locations under my home directory
#-----------------------------------------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

#-----------------------------------------------------------------------------------------
# Folder locations under my working directory
#-----------------------------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
loc.MRPRESSO <- paste0(workingDir,"MR_ICC_GSCAN_201806/MR-PRESSO/")
loc.MRPRESSO.output <- paste0(loc.MRPRESSO,"output")
loc.MRPRESSO.tabulated <- paste0(loc.MRPRESSO,"result-tabulated/")
#dir.create(loc.MRPRESSO.tabulated)

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

# A list of file paths to subset
list.files(path=loc.MRPRESSO.output
           ,full.names = TRUE)

# Get a vector of file paths of all MR result files in the input folder
pattern.file.names <- "MRPRESSO-global-test.*\\.tsv"

file.paths <- list.files(path=loc.MRPRESSO.output
                         ,pattern = pattern.file.names
                         ,full.names = TRUE) # length(file.paths) 4

#--------------------------------------------------------------------------------------------------
# Combine all MRPRESSO global test files
#--------------------------------------------------------------------------------------------------
# Scale up causal estimate of caffeine exposure up using 1 cup change of coffee (75mg) or tea (40 m) consumption
caffeine.per.cup.coffee <- 75
caffeine.per.cup.tea <- 40
caffeine.per.cup <- (caffeine.per.cup.coffee+caffeine.per.cup.tea)/2

base.MRPRESSO.global.outlier.corrected <- data.frame()

for (i in 1:length(file.paths)){
  filePath <- file.paths[i]
  # Import MR-PRESSO global test result file
  # Convert effect size estimates and 95% confidence intervals
  MRPRESSO.global.outlier.corrected <- ImportATabSeparatedFile(input.file.path = filePath
                                                               , data.name = "MRPRESSO.global") %>% 
    tidyr::separate(col = "exposure", into=c("exposure.consortium","exposure"),sep="-",remove=TRUE) %>%
    tidyr::separate(col = "outcome", into=c("outcome.consortium","outcome"),sep="-",remove=TRUE) %>%
    dplyr::rename( exposure.clumping.p1 = clumping.p1
                  ,pval=p.value) %>%
    dplyr::mutate(method="MRPRESSO"
                  ,var.type.exposure=dplyr::case_when( exposure %in% c("SC","SI","CI") ~ "binary"
                                           ,exposure %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                  ,var.type.outcome=dplyr::case_when( outcome %in% c("SC","SI","CI") ~ "binary"
                                           ,outcome %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                  ,b=causal.estimate
                  ,effect.size.estimate=dplyr::case_when( 
                     var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b
                    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ b*caffeine.per.cup
                    ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b)
                    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(b*caffeine.per.cup)
                    ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*b
                    ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*b)
                                                  )
                  ,effect.size.LBound=dplyr::case_when(
                     var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b-1.96*sd
                    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b-1.96*sd)
                    ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b-1.96*sd)
                    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b-1.96*sd))
                    ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b-1.96*sd)
                    ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b-1.96*sd))
                                                )
                  ,effect.size.UBound=dplyr::case_when(
                     var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b+1.96*sd
                    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b+1.96*sd)
                    ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b+1.96*sd)
                    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b+1.96*sd))
                    ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b+1.96*sd)
                    ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b+1.96*sd))
                                                )
                  ) %>%
    dplyr::select_(.dots=c("MR.analysis","exposure.consortium","exposure","exposure.clumping.p1","outcome.consortium","outcome","method","b","sd","pval","effect.size.estimate","effect.size.LBound","effect.size.UBound"))
  
  # Append result to the base data set
  base.MRPRESSO.global.outlier.corrected <- rbind( base.MRPRESSO.global.outlier.corrected
                                                  ,MRPRESSO.global.outlier.corrected)
  
}

# dim(base.MRPRESSO.global.outlier.corrected) 8 13

ExportFileTabSeparated(data = base.MRPRESSO.global.outlier.corrected
                       ,missing.values.as = ""
                       , output.file.path = paste0(loc.MRPRESSO.tabulated,"MR-PRESSO-global-test-results",".tsv"))






#--------------------------------------------------------------------------------------------------
# Calculate odds ratios and 95% CI using coefficients from MR-Egger, weighted median, weighted mode
#--------------------------------------------------------------------------------------------------
tem1 <- ImportATabSeparatedFile(input.file.path = "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/two-sample-MR/result-tabulated/MR-analysis-results.tsv"
                        ,data.name = "MR.analysis.results") %>% 
  dplyr::filter(exposure.consortium=="UKB" & outcome.consortium=="ICC" & outcome=="CI" & method.abbreviated %in% c("Egger","W Median","W Mode")) %>%
  dplyr::filter((exposure=="CCPD" & exposure.clumping.p1==1e-05) | (exposure=="CCPD" & exposure.clumping.p1==5e-08) | (exposure=="ESDPW" & exposure.clumping.p1==1e-05) | (exposure=="PYOS" & exposure.clumping.p1==1e-05)) %>%
  dplyr::mutate(odds.ratio=exp(b)
                ,odds.ratio.lower.bound=exp(b - 1.96*se)
                ,odds.ratio.upper.bound=exp(b + 1.96*se)) %>%
  dplyr::select_(.dots=colnames(base.MRPRESSO.global.outlier.corrected)) # dim(tem1) 12 10

#----------------------------------------------------------------------------
# Combine the two data.frames
#----------------------------------------------------------------------------
tem2 <- rbind(base.MRPRESSO.global.outlier.corrected,tem1) # dim(tem2) 16 10

# Export result file. Set missing values to white space for SAS
ExportFileTabSeparated(data = tem2
                       ,missing.values.as = ""
                       , output.file.path = paste0(loc.MRPRESSO.tabulated,"odds-ratio_exposure-UKB-CCPD-ESDPW-PYOS_outcome-ICC-CI_MR-sensitivity-analyses",".tsv"))

#setwd(locScripts)
#file.copy("MR_step06-05-05_tabulate-MRPRESSO-results.R")
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#