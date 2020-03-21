#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-07_run_heterogeneity-test.R
# Modified from : MR_step06-03-06_run_MR-leave-one-out_MR-IVW.R
# Date created  : 20190822
# Purpose       : Calculate heterogeneity statistics on harmonised data using mr_heterogeneity()
#                 Merge two-sample MR result with heterogeneity test result as a single file
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency:     
# Function external:  
# Type  Files
#----------------------------------------------------------------------------------------
# Input paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# Input Sys.glob(path=paste0(loc.twoSampleMR.harmonised,"harmonised-data*\\.tsv")) # 220 files

# Outpu paste0(loc.twoSampleMR.tabulated,"heterogeneity-test_two-sample-MR-results_sorted.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190822  Exported paste0(loc.twoSampleMR.tabulated,"heterogeneity-test_two-sample-MR-results_sorted.tsv")
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
loc.twoSampleMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/two-sample-MR/")
loc.twoSampleMR.harmonised <- paste0(loc.twoSampleMR,"input/harmonised-data/")
loc.twoSampleMR.tabulated <- paste0(loc.twoSampleMR,"result-tabulated/")
# Create output folder
loc.twoSampleMR.heterogeneity <- paste0(loc.twoSampleMR,"output_heterogeneity-test/")
#dir.create(loc.twoSampleMR.heterogeneity)

#-----------------------------------------------------------------------------------------
# Load packages, import external functions
#-----------------------------------------------------------------------------------------
library("TwoSampleMR",lib.loc="/software/R/R-3.5.1/lib64/R/library")
library("dplyr",lib.loc="/software/R/R-3.4.1/lib64/R/library")
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#-----------------------------------------------------------------------------------------
# Import data
#-----------------------------------------------------------------------------------------
# ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
#                         ,data.name = "two.sample.MR.results") # dim(two.sample.MR.results) 379 21

#---------------------------------------------------
# Get a list of paths of harmonised data files
#---------------------------------------------------
# Method 1: use list.files() and glob2rx() 
## glob2rx() uses * as a wildcard and expands the string into regular expression for R. Note regex in R is different from regex in Linux
## Here glob2rx() added ^, dot, \\ and $
### regex in R: (1) dot usually means "match any character" (2) * means The preceding item will be matched zero or more times. .* is the same as * in linux
pattern.file.names <- glob2rx("harmonised-data_exposure*.tsv") # "^harmonised-data_exposure.*\\.tsv$"

file.paths.harmonised.data <- list.files(path=loc.twoSampleMR.harmonised
                                         ,pattern = pattern.file.names
                                         ,full.names = TRUE) # length(file.paths.harmonised.data) 220
# Method 2: use Sys.glob()
# Sys.glob() Expands wildcard on file paths
file.paths.harmonised.data <- Sys.glob(path=paste0(loc.twoSampleMR.harmonised,"harmonised-data*\\.tsv")) # length(file.paths.harmonised.data) 220

#---------------------------------------------------
# Run heterogeneity test on all harmonised data
#---------------------------------------------------
base.heterogeneity <- data.frame()

# Loop thru each file path, calculating heterogeneity statistics 
## Number of iterations: 220
for (i in 1:length(file.paths.harmonised.data)){
  file.path <- file.paths.harmonised.data[i]
  ImportATabSeparatedFile(input.file.path = file.path, data.name = "data.harmonised") # dim(data.harmonised) 63 28
  #----------------------------------------------------------------------------------------------
  # Run heterogeneity test
  #----------------------------------------------------------------------------------------------
  heterogeneity <- TwoSampleMR::mr_heterogeneity(data.harmonised
                                                 ,method_list=c("mr_egger_regression", "mr_ivw"))
  # Create output file name by changing the name prefix to the pattern
  output.file.name <- gsub(basename(file.path)
                           ,pattern = "harmonised-data"
                           ,replacement = "heterogeneity-test")
  
  # Export test result as a single file
  ExportFileTabSeparated(data=heterogeneity
                         ,missing.values.as = "NA"
                         , output.file.path = paste0(loc.twoSampleMR.heterogeneity,output.file.name))
  
  base.heterogeneity <- rbind(heterogeneity,base.heterogeneity)

} # dim(base.heterogeneity) 410 8

# # Change factor columns to character for merging
# base.heterogeneity[,c("exposure","outcome","method")] <- lapply(base.heterogeneity[,c("exposure","outcome","method")]
#                                                                 ,as.character)
# 
# # Horizontally combine two sample MR results (left table) and heterogeneity test result
# hori.merge <- dplyr::left_join(two.sample.MR.results
#                                ,base.heterogeneity
#                                ,by=c(  "id.exposure"="id.exposure"
#                                       ,"exposure"="exposure"
#                                      ,"id.outcome"="id.outcome"
#                                      ,"outcome"="outcome"
#                                      ,"method"="method")) # dim(hori.merge) 657 24
# 
# # Sort data by this order 
# ## Ascending exposure.consortium, ascending exposure, ascending outcome_consortium, ascending outcome, DESCENDING exposure.clumping.p1 and ## MR methods: IVW, MR Egger, W Mode, W median
# hori.merge$method.abbreviated <- factor(hori.merge$method.abbreviated, level=c("IVW","Egger","W Median","W Mode"))
# 
# hori.merge.sorted <- hori.merge %>% 
#   dplyr::arrange(exposure.consortium
#                  ,exposure
#                  ,outcome.consortium
#                  ,outcome
#                  ,desc(exposure.clumping.p1)
#                  ,method.abbreviated) %>%
#   dplyr::mutate(row.order=c(1:nrow(.)))
# 
# 
# # Export data of all test results to a single file
# ExportFileTabSeparated(data = hori.merge.sorted
#                        ,missing.values.as = "" # set missing to blank for SAS
#                        ,output.file.path = paste0(loc.twoSampleMR.tabulated,"heterogeneity-test_two-sample-MR-results_sorted.tsv"))
# 

setwd(locScripts)
file.copy( from="MR_step06-03-07_run_heterogeneity-test.R"
          ,to="MR_step06-03-08_combine-results_two-sample-MR_heterogeneity-test.R")

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#


#------------------------------------------------------------------------------
# Remove 2 SNPs (rs2472297, rs4410790) and then run inverse-variance weighted MR
#------------------------------------------------------------------------------
# Exclude 2 SNPs (rs2472297, rs4410790)
harmonised.UKBCCPD.ICCCI.5e8.rm2SNPs <- subset(harmonised.UKBCCPD.ICCCI.5e8, !SNP %in% c("rs2472297","rs4410790")) # dim(harmonised.UKBCCPD.ICCCI.5e8.rm2SNPs) 15 28

# Perform IVW
IVW.harmonised.UKBCCPD.ICCCI.5e8.rm2SNPs <- TwoSampleMR::mr(harmonised.UKBCCPD.ICCCI.5e8.rm2SNPs
                                                            , method_list=c("mr_ivw"))

# Analysis result:
# > IVW.harmonised.UKBCCPD.ICCCI.5e8.rm2SNPs
# id.exposure id.outcome outcome exposure                    method nsnp           b        se      pval
# 1      s6xdUH     1HvR1n      CI     CCPD Inverse variance weighted   15 -0.05617675 0.0926269 0.5441933

