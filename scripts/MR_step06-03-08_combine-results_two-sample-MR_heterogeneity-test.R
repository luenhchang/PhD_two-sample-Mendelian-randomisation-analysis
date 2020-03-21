#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-07_run_heterogeneity-test.R
# Modified from : MR_step06-03-06_run_MR-leave-one-out_MR-IVW.R
# Date created  : 20190926
# Purpose       : Merge two-sample MR result with heterogeneity test result as a single file
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency: MR_step06-03-05_tabulate-two-sample-MR-analysis-results.R    
# Function external:  
# Type  Files
#----------------------------------------------------------------------------------------
# Input paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# Input Sys.glob(path=paste0(loc.twoSampleMR.harmonised,"harmonised-data*\\.tsv")) # 220 files

# Outpu paste0(loc.twoSampleMR.tabulated,"heterogeneity-test_two-sample-MR-results_sorted.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20191024, 20191010, 20191002, 20190927, 20190822  
# Exported heterogeneity-test_two-sample-MR-results_sorted.tsv
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
#loc.twoSampleMR.harmonised <- paste0(loc.twoSampleMR,"input/harmonised-data/")
loc.twoSampleMR.tabulated <- paste0(loc.twoSampleMR,"result-tabulated/")
loc.twoSampleMR.heterogeneity <- paste0(loc.twoSampleMR,"output_heterogeneity-test")
#dir.create(loc.twoSampleMR.heterogeneity)

#-----------------------------------------------------------------------------------------
# Load packages, import external functions
#-----------------------------------------------------------------------------------------
# library("TwoSampleMR",lib.loc="/software/R/R-3.5.1/lib64/R/library")
# library("dplyr",lib.loc="/software/R/R-3.4.1/lib64/R/library")
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#-----------------------------------------------------------------------------------------
# Import data
#-----------------------------------------------------------------------------------------
ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
                        ,data.name = "two.sample.MR.results") # dim(two.sample.MR.results) 657 22

#---------------------------------------------------
# Combine heterogeneity test result files as a single data
#---------------------------------------------------
pattern.file.names <- glob2rx("heterogeneity-test*.tsv") # "^harmonised-data_exposure.*\\.tsv$"

file.paths.heterogeneity.test <- list.files(path=loc.twoSampleMR.heterogeneity
                                         ,pattern = pattern.file.names
                                         ,full.names = TRUE) # length(file.paths.heterogeneity.test) 220

base.heterogeneity <- data.frame()

# Loop thru each file path, calculating heterogeneity statistics 
## Number of iterations: 220
for (i in 1:length(file.paths.heterogeneity.test)){
  file.path <- file.paths.heterogeneity.test[i]
  tryCatch({

  # Read herterogeneity test result tsv files.
  ## 10 of 220 files have 1 line, the other 210 files have 2 lines of data (not including header rows)  
    ImportATabSeparatedFile(input.file.path =file.path, data.name = "data.heterogeneity") 
    base.heterogeneity <- rbind(base.heterogeneity,data.heterogeneity)
    # End tryCatch() function
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
} # dim(base.heterogeneity) 410 8


# Change factor columns to character for merging
base.heterogeneity[,c("exposure","outcome","method")] <- lapply(base.heterogeneity[,c("exposure","outcome","method")]
                                                                ,as.character)

# Horizontally combine two sample MR results (left table) and heterogeneity test result
hori.merge <- dplyr::left_join(two.sample.MR.results
                               ,base.heterogeneity
                               ,by=c(  "id.exposure"="id.exposure"
                                       ,"exposure"="exposure"
                                       ,"id.outcome"="id.outcome"
                                       ,"outcome"="outcome"
                                       ,"method"="method")) # dim(hori.merge) 657 25

# Sort data by this order 
## Ascending exposure.consortium, ascending exposure, ascending outcome_consortium, ascending outcome, DESCENDING exposure.clumping.p1 and ## MR methods: IVW, MR Egger, W Mode, W median
hori.merge$method.abbreviated <- factor(hori.merge$method.abbreviated, level=c("IVW","Egger","W Median","W Mode"))

hori.merge.sorted <- hori.merge %>% 
  dplyr::arrange(exposure.consortium
                 ,exposure
                 ,outcome.consortium
                 ,outcome
                 ,desc(exposure.clumping.p1)
                 ,method.abbreviated) %>%
  dplyr::mutate(row.order=c(1:nrow(.))) # dim(hori.merge.sorted) 657 26

# Export data of all test results to a single file
ExportFileTabSeparated(data = hori.merge.sorted
                       ,missing.values.as = "" # set missing to blank for SAS
                       ,output.file.path = paste0(loc.twoSampleMR.tabulated,"heterogeneity-test_two-sample-MR-results_sorted.tsv"))

setwd(locScripts)

# Multiple testing threshold is calculated as 2 (unique exposure sample and substance) * 6 (unique outcome sample and substance) = 12. Manually copy this threshold to NU4tabSup04_heterogeneity-test_two-sample-MR-results.sas
hori.merge.sorted.want <- hori.merge.sorted %>% dplyr::filter(exposure.outcome.status=="reasonable two-sample MR") # dim(hori.merge.sorted.want) 164 26
unique.exposure.sample.substance.groups <- unique(hori.merge.sorted.want[,c("exposure.consortium","exposure.substance")]) # nrow(unique.exposure.sample.substance.groups) 2

unique.outcome.sample.substance.groups <- unique(hori.merge.sorted.want[,c("outcome.consortium","outcome.substance")]) # nrow(unique.outcome.sample.substance.groups) 6 

# file.copy( from="MR_step06-03-07_run_heterogeneity-test.R"
#           ,to="MR_step06-03-08_combine-results_two-sample-MR_heterogeneity-test.R")

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

