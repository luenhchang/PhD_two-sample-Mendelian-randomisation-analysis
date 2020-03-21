##################################################################################
# filename: MR_step06-05-01_prepare-input-files-for-MR-PRESSO.R
# program author: Chang
# purpose: prepare input files and a tsv file to loop through for running MR-PRESSO at next step
# date created: 20190329
# Function internal: 
# Function external: ImportATabSeparatedFile()
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Input paste0(loc.harmonised.data,"/harmonised-data_exposure-clumped-dpw-noICC-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5-linear-BETA-added_outcome-GSCAN-CPD.tsv")
# Input paste0(loc.harmonised.data,"/harmonised-data_exposure-clumped-GWAS-UKB-caffeine-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5_outcome-UKB-CPD.tsv")
# Input paste0(loc.harmonised.data,"/harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5_outcome-UKB-caffeine.tsv")
# Input paste0(loc.harmonised.data,"/harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine.tsv")

# Ouput paste0(loc.MRPRESSO.input,"/harmonised-data_exposure-clumped-dpw-noICC-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5-linear-BETA-added_outcome-GSCAN-CPD.tsv")
# Ouput paste0(loc.MRPRESSO.input,"/harmonised-data_exposure-clumped-GWAS-UKB-caffeine-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5_outcome-UKB-CPD.tsv")
# Ouput paste0(loc.MRPRESSO.input,"/harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5_outcome-UKB-caffeine.tsv")
# Ouput paste0(loc.MRPRESSO.input,"/harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine.tsv")

# Ouput paste0(loc.MRPRESSO.input,"file-info_harmonised-data-selected-pairs-exposures-outcomes",".tsv")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20191011  Exported file-info_harmonised-data-selected-pairs-exposures-outcomes.tsv
# 20190822  Added one trait pair. Exported file-info_harmonised-data-selected-pairs-exposures-outcomes.tsv
# 20190818  Exported file-info_harmonised-data-selected-pairs-exposures-outcomes.tsv
# 20190714  Exported file-info_harmonised-data-selected-pairs-exposures-outcomes.tsv
# 20190419  Exported the 4 files above
# 20190329  Exported the 4 files above
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Folder locations under my home directory
#-----------------------------------------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

#-----------------------------------------------------------------------------------------
# Folder locations under my working directory
#-----------------------------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
loc.MR <- paste0(workingDir,"MR_ICC_GSCAN_201806/")
loc.harmonised.data <- paste0(loc.MR,"two-sample-MR/input/harmonised-data")
loc.MRPRESSO <- paste0(loc.MR,"MR-PRESSO/")
loc.MRPRESSO.input <- paste0(loc.MRPRESSO,"input/")
loc.MRPRESSO.input.archive <- paste0(loc.MRPRESSO.input,"archive")
loc.MRPRESSO.output <- paste0(loc.MRPRESSO,"output/")
# dir.create(loc.MRPRESSO)
# dir.create(loc.MRPRESSO.input)
# dir.create(loc.MRPRESSO.input.archive)
# dir.create(loc.MRPRESSO.output)

#------------------------------------------------------------------
# Install package
#------------------------------------------------------------------
# https://www.nature.com/articles/s41588-018-0099-7

if (!require("devtools")) { install.packages("devtools") } else {}

devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO,lib.loc = "/software/R/R-3.5.1/lib64/R/library")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#-----------------------------------------------------------------
# Copy selected harmonised data files to the input folder
#-----------------------------------------------------------------

# Exposure                          clumping p1 Outcome
#-------------------------------------------------------------------------------
# GSCAN regular smoking initiation  p1-1e-5     UKB-caffeine
# GSCAN regular smoking initiation  p1-5e-8     UKB-caffeine
# ICC cannabis initiation           p1-5e-8     GSCAN regular smoking initiation
# ICC cannabis initiation           p1-5e-8     UKB pack years of smoking
#-------------------------------------------------------------------------------

# Specify patterns in file names to search in file paths. Note file paths do not begin with file name. 
# pattern.1 <- "harmonised-data_exposure-clumped-dpw-noICC-LDWindow-kb-10000-R2-0\\.01-p1-1e-5-p2-1e-5-linear-BETA-added_outcome-GSCAN-CPD\\.tsv$"
# pattern.2 <- "harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0\\.01-p1-1e-5-p2-1e-5_outcome-UKB-caffeine\\.tsv$"
# pattern.3 <- "harmonised-data_exposure-clumped-GWAS-UKB-caffeine-LDWindow-kb-10000-R2-0\\.01-p1-1e-5-p2-1e-5_outcome-UKB-CPD\\.tsv$"
# pattern.4 <- "harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0\\.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine\\.tsv$"
# pattern.5 <- "harmonised-data_exposure-clumped-GWAS-UKB-caffeine-LDWindow-kb-10000-R2-0\\.01-p1-5e-8-p2-1e-6_outcome-ICC-CI\\.tsv$"

pattern.1 <- "harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0\\.01-p1-1e-5-p2-1e-5_outcome-UKB-caffeine\\.tsv$"
pattern.2 <- "harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0\\.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine\\.tsv$"
pattern.3 <- "harmonised-data_exposure-clumped-GWAS-ICC-CI-LDWindow-kb-10000-R2-0\\.01-p1-5e-8-p2-1e-6_outcome-GSCAN-SI\\.tsv$"
pattern.4 <- "harmonised-data_exposure-clumped-GWAS-ICC-CI-LDWindow-kb-10000-R2-0\\.01-p1-5e-8-p2-1e-6_outcome-UKB-PYOS\\.tsv$"

# Combine multiple patterns as one string
patterns <- paste0(pattern.1,"|",pattern.2,"|",pattern.3,"|",pattern.4)

# Get paths of the files method 1
source.file.paths <- list.files(path = loc.harmonised.data
                                ,pattern = patterns
                                ,full.names = TRUE) # length(source.file.paths) 4
# Get paths of the files method 2
source.file.paths <- grep(patterns
                        ,list.files(path= loc.harmonised.data
                                    ,full.names = T)
                        ,value = TRUE) # length(source.file.paths) 4

# Create destination file paths
destin.file.paths <- paste0(loc.MRPRESSO.input,basename(source.file.paths)) # length(destin.file.paths) 4

# Manually move existing files to the archive folder
file.copy(from= source.file.paths
          ,to= destin.file.paths
          ,overwrite = TRUE)

#--------------------------------------------------------------------------------
# Create a file for looping through each row for running MR-PRESSO at next step
## MR-PRESSO global test evaluates for the presence of horizontal pleiotropy
## MR-PRESSO outlier test identifies outlier SNP, resulting in 1 row per SNP data
#--------------------------------------------------------------------------------
# harmonised.data.file.paths <- c(paste0(loc.MRPRESSO.input,file.CCPD.1e5.CI)
#                                 ,paste0(loc.MRPRESSO.input,file.CCPD.5e8.CI)
#                                 ,paste0(loc.MRPRESSO.input,file.ESDPW.1e5.CI)
#                                 ,paste0(loc.MRPRESSO.input,file.PYOS.1e5.CI))

#numb.harmonised.data.files <- length(harmonised.data.file.paths)
numb.harmonised.data.files <- length(destin.file.paths)

signi.thres.Bonferroni <- 0.05/numb.harmonised.data.files

# Manually list the trait pairs as the order of the destination file paths. Be aware that data have been sorted
# Exposure              clumping p1 Outcome
#------------------------------------------------
# GSCAN drinks per week p1-1e-5     GSCAN CPD
# UKB-caffeine          p1-1e-5     UKB-CPD
# UKB-caffeine          p1-5e-8     ICC-CI
# GSCAN SI              p1-1e-5     UKB-caffeine
# GSCAN SI              p1-5e-8     UKB-caffeine
#------------------------------------------------

files.to.run.MRPRESSO <- data.frame(file.path.harmonised=destin.file.paths
                                    ,exposure= c("ICC-CI","ICC-CI","GSCAN-SI","GSCAN-SI")
                                    ,clumping.p1=c("5e-8","5e-8","1e-5","5e-8")
                                    ,outcome= c("GSCAN-SI","UKB-PYOS","UKB-caffeine","UKB-caffeine") 
                                    ,colname_for_SNP=rep("SNP", numb.harmonised.data.files)
                                    ,colname_for_beta_outcome=rep("beta.outcome",numb.harmonised.data.files)
                                    ,colname_for_beta_exposure=rep("beta.exposure",numb.harmonised.data.files)
                                    ,colname_SD_outcome=rep("se.outcome",numb.harmonised.data.files)
                                    ,colname_SD_exposure=rep("se.exposure",numb.harmonised.data.files)
                                    ,numb_simulation=rep("10000",numb.harmonised.data.files)
                                    ,significance_threshold=rep(signi.thres.Bonferroni,numb.harmonised.data.files)
                                    ,stringsAsFactors = F) # dim(files.to.run.MRPRESSO) 4 11

ExportFileTabSeparated( data = files.to.run.MRPRESSO
                        ,output.file.path = paste0(loc.MRPRESSO.input
                                                   ,"file-info_harmonised-data-selected-pairs-exposures-outcomes"
                                                   ,".tsv"))

#file.copy( paste0(locScripts,"MR_step06-05-01_prepare-input-files-for-MR-PRESSO.R"),paste0(locScripts,"MR_step06-05-02_job-script_run-MR-PRESSO.R"))

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#
