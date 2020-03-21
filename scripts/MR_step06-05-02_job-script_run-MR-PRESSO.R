#!/usr/bin/env Rscript

##################################################################################
# Filename: MR_step06-05-02_job-script_run-MR-PRESSO.R
# Modified from: MR_step06-05-01_prepare-input-files-for-MR-PRESSO.R
# program author: Chang
# purpose: Being submitted as a PBS job for running MR-PRESSO
# date created: 20190419
# Function internal: 
# Function external: 
# Note: every argument is passed as a string. Convert variable to numeric if the function works on numeric value
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190419  
#-----------------------------------------------------------------------------------------

homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
loc.MR <- paste0(workingDir,"MR_ICC_GSCAN_201806/")
loc.harmonised.data <- paste0(loc.MR,"two-sample-MR/input/harmonised-data/")
loc.MRPRESSO <- paste0(loc.MR,"MR-PRESSO/")
loc.MRPRESSO.input <- paste0(loc.MRPRESSO,"input/")
loc.MRPRESSO.output <- paste0(loc.MRPRESSO,"output/")

# Get arguments specified in part 2
arguments.passed.bash.to.R <- commandArgs(trailingOnly = TRUE)
print(paste0(arguments.passed.bash.to.R))

# Check if the arguments contain nothing 
if (length(arguments.passed.bash.to.R) < 1)
  stop("Missing argument: num.rows")

#--------------------------------------------------------
# Pass command line arguments to variables in R 
#--------------------------------------------------------
print("==================Passing command line arguments sequentially to variables in R==================")
R.variables <- c("file.path.harmonised.data"
                 ,"exposure.trait"
                 ,"outcome.trait"
                 ,"colname.SNP"
                 ,"colname.beta.outcome"
                 ,"colname.beta.exposure"
                 ,"colname.SD.outcome"
                 ,"colname.SD.exposure"
                 ,"numb.simulation"
                 ,"significance.threshold"
                 ,"clumping.p1"
                 ,"output.folder.path") # length(R.variables) 12

# Pass arguments to variables. The length of R.variables and the argument list must be the same
for (i in 1:length(R.variables)){
  print(paste0(R.variables[i]," value=",arguments.passed.bash.to.R[[i]]))
  assign(R.variables[i],arguments.passed.bash.to.R[[i]])
}

#----------------------------------------------------------    
# Test code with fixed values for the variables
## Comment this part out when the code is working fully
#----------------------------------------------------------
# file.path.harmonised.data <- "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/MR-PRESSO/input/harmonised-data_exposure-clumped-GWAS-UKB-CCPD-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-ICC-CI"
# exposure.trait <- "UKB-CCPD"
# outcome.trait <- "ICC-CI"
# colname.SNP <- "SNP"
# colname.beta.outcome <- "beta.outcome"
# colname.beta.exposure <- "beta.exposure"
# colname.SD.outcome <- "se.outcome"
# colname.SD.exposure <- "se.exposure"
# numb.simulation <- 10000
# significance.threshold <- 0.05
# clumping.p1 <- "5e-8"
# output.folder.path <- "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/MR-PRESSO/output/"

#------------------------------------------------------------------
# Install package
# Import functions
#------------------------------------------------------------------
# if (!require("devtools")) { install.packages("devtools") } else {}
# 
# devtools::install_github("rondolab/MR-PRESSO")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#-----------------------------------------------------------------
# Import a harmonised data 
#-----------------------------------------------------------------
ImportATabSeparatedFile(input.file.path = file.path.harmonised.data
                        ,data.name = "harmonised.data") # dim(harmonised.CCPD.5e8.CI) 17 28

# Sort input data by SNP Rs ID
data <- harmonised.data
data.sorted <- data[order(data[,colname.SNP]),]

#--------------------------------------------------------------------------------
# Run Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO) test
## MR-PRESSO global test evaluates for the presence of horizontal pleiotropy
## MR-PRESSO outlier test identifies outlier SNP, resulting in 1 row per SNP data
#--------------------------------------------------------------------------------
MRPRESSO.results <- MRPRESSO::mr_presso(data=data.sorted
                                          , BetaOutcome = colname.beta.outcome
                                          , BetaExposure = colname.beta.exposure
                                          , SdOutcome = colname.SD.outcome
                                          , SdExposure = colname.SD.exposure
                                          , OUTLIERtest = TRUE
                                          , DISTORTIONtest = TRUE
                                          , NbDistribution = as.numeric(numb.simulation)
                                          , SignifThreshold = as.numeric(significance.threshold))
  
#---------------------------------------------------------------------------
# ---------------- Tabulate Get MR-PRESSO global test result
#---------------------------------------------------------------------------
MRPRESSO.global <- MRPRESSO.results[[1]] # dim(MRPRESSO.global) 2 6

# Rename columns. There are space in some column names
colnames(MRPRESSO.global) <- c("exposure","MR.analysis","causal.estimate","sd","T.stat","p.value")

# Replace default exposure, add outcome name
MRPRESSO.global$exposure <- exposure.trait
MRPRESSO.global$clumping.p1 <- clumping.p1
MRPRESSO.global$outcome <- outcome.trait

# Export MR-PRESSO global test result
file.prefix <- "MRPRESSO-global-test_"
file.suffix.exposure.outcome <- paste0("exposure-",exposure.trait,"-",clumping.p1,"_outcome-",outcome.trait)
ExportFileTabSeparated( data = MRPRESSO.global
                        ,output.file.path = paste0(output.folder.path
                                                   ,file.prefix
                                                   ,file.suffix.exposure.outcome
                                                   ,".tsv"))
#---------------------------------------------------------------------------
# ---------------- Tabulate outlier test result
#---------------------------------------------------------------------------
## data dimension: 1 SNP per row
MRPRESSO.outliers <- data.frame(MRPRESSO.results[[2]]["Outlier Test"]
                                ,stringsAsFactors = F) # dim(MRPRESSO.outliers) 176 2
# Merge SNP Rs ID and outlier test results
MRPRESSO.outliers2 <- cbind(data.frame(SNP=data.sorted[,colname.SNP])
                            ,MRPRESSO.outliers)

# Add exposure, clumping p1, outcome name
MRPRESSO.outliers2$Exposure <- exposure.trait
MRPRESSO.outliers2$clumping.p1 <- clumping.p1
MRPRESSO.outliers2$outcome <- outcome.trait

# Export MR-PRESSO outlier test result
file.prefix <- "MRPRESSO-outlier-test_"
file.suffix.exposure.outcome <- paste0("exposure-",exposure.trait,"-",clumping.p1,"_outcome-",outcome.trait)
ExportFileTabSeparated( data = MRPRESSO.outliers2
                        ,output.file.path = paste0(output.folder.path
                                                   ,file.prefix
                                                   ,file.suffix.exposure.outcome
                                                   ,".tsv"))

# MR-PRESSO Outlier indices (which SNP is the outlier?)
outliers.indices <- unlist(MRPRESSO.results[[2]]["Distortion Test"])[1]
outliers.SNPs <- data.sorted[,colname.SNP][outliers.indices]

MRPRESSO.outlier.indices <- data.frame(exposure.trait=exposure.trait
                                       ,clumping.p1=clumping.p1
                                       ,outcome.trait=outcome.trait
                                       ,distortion.test.outlier.indices=outliers.indices
                                       ,distortion.test.outlier.SNP=outliers.SNPs
                                       ,stringsAsFactors = F)

# Export MR-PRESSO outlier test result
file.prefix <- "MRPRESSO-distortion-test_"
file.suffix.exposure.outcome <- paste0("exposure-",exposure.trait,"-",clumping.p1,"_outcome-",outcome.trait)
ExportFileTabSeparated( data = MRPRESSO.outlier.indices
                        ,output.file.path = paste0(output.folder.path
                                                   ,file.prefix
                                                   ,file.suffix.exposure.outcome
                                                   ,".tsv"))

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#
