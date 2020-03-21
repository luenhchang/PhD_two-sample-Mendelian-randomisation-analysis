#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step06-03-02_jobScript-harmonise-exposure-outcome.R
# Modified from : PRS_UKB_201711_step21-05-01_jobScript_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.R
# Date created  : 20190405
# Purpose       : Harmonise exposure (clumped GWAS) and outcomes (QCed GWASs) using twoSampleMR package in a R script (the current script) through qsub (see MR_step06-03-03) a bash script (see MR_step06-03-04)
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency:     
# Function external:  Run.two.sample.Mendelian.randomisation()
# Type  Files
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# Locations of input, output folders
#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

print("============================================================================================")
print("=================== Running the R script from this line belows =============================")
print("============================================================================================")

# Get arguments, specified in PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh, from commandline to R 
arguments.passed.from.bash.to.R <- commandArgs(trailingOnly = TRUE)

# Check if the arguments contain nothing 
if (length(arguments.passed.from.bash.to.R) < 1)
  stop("Missing argument: num.rows")

# Extract individual elements of the argument 
print("============================================================================================")
print("======================Passing arguments from Shell to R===============================")
print("============================================================================================")

# ArguNum Corresponding Shell variable at MR_step06-03-03
#-----------------------------------------------------
# 1       ${tsv_1_filePath} 
# 2       ${tsv_1_consortium} 
# 3       ${tsv_1_trait} 
# 4       ${tsv_1_file_delimiter}
# 5       ${tsv_1_colname_for_SNP} 
# 6       ${tsv_1_colname_for_beta} 
# 7       ${tsv_1_colname_for_SE} 
# 8       ${tsv_1_colname_for_effect_allele} 
# 9       ${tsv_1_colname_for_other_allele} 
# 10      ${tsv_1_colname_for_P_value} 
# 11      ${tsv_1_clumping_p1_value} 
# 12      ${tsv_1_clumping_p2_value}
# 13      ${tsv_2_filePath}
# 14      ${tsv_2_consortium} 
# 15      ${tsv_2_trait} 
# 16      ${tsv_2_file_delimiter}
# 17      ${tsv_2_colname_for_SNP} 
# 18      ${tsv_2_colname_for_beta} 
# 19      ${tsv_2_colname_for_SE} 
# 20      ${tsv_2_colname_for_effect_allele} 
# 21      ${tsv_2_colname_for_other_allele} 
# 22      ${tsv_2_colname_for_P_value} 
# 23      ${iteration}
# 24      ${filePath_output_harmonised_data} 
# 25      ${filePath_output_MR_analysis}
# 26      ${folderPath_output_MR_report}
# 27      ${filePath_output_hori_pleio}
#-----------------------------------------------------

# Pass Shell variables related to exposure to R
here.exposure.file.path <- arguments.passed.from.bash.to.R[[1]]
here.exposure.consortium <- arguments.passed.from.bash.to.R[[2]]
here.exposure.trait <- arguments.passed.from.bash.to.R[[3]]
here.exposure.file.delimiter <- arguments.passed.from.bash.to.R[[4]]

here.exposure.file.colname.for.SNP <- arguments.passed.from.bash.to.R[[5]]  
here.exposure.file.colname.for.beta <- arguments.passed.from.bash.to.R[[6]]  
here.exposure.file.colname.for.standard.error <- arguments.passed.from.bash.to.R[[7]]
here.exposure.file.colname.for.effect.allele <- arguments.passed.from.bash.to.R[[8]]
here.exposure.file.colname.for.other.allele <- arguments.passed.from.bash.to.R[[9]]
here.exposure.file.colname.for.P.value <- arguments.passed.from.bash.to.R[[10]]
here.exposure.clumping.p1 <- arguments.passed.from.bash.to.R[[11]]
here.exposure.clumping.p2 <- arguments.passed.from.bash.to.R[[12]]

# Pass Shell variables related to outcome to R
here.outcome.file.path <- arguments.passed.from.bash.to.R[[13]]
here.outcome.consortium <- arguments.passed.from.bash.to.R[[14]]
here.outcome.trait <- arguments.passed.from.bash.to.R[[15]]
here.outcome.file.delimiter <- arguments.passed.from.bash.to.R[[16]]
here.outcome.file.colname.for.SNP <- arguments.passed.from.bash.to.R[[17]]
here.outcome.file.colname.for.beta <- arguments.passed.from.bash.to.R[[18]]
here.outcome.file.colname.for.standard.error <- arguments.passed.from.bash.to.R[[19]]
here.outcome.file.colname.for.effect.allele <- arguments.passed.from.bash.to.R[[20]]
here.outcome.file.colname.for.other.allele <- arguments.passed.from.bash.to.R[[21]]
here.outcome.file.colname.for.P.value <- arguments.passed.from.bash.to.R[[22]]

# Pass Shell variables related to others
here.iteration <- arguments.passed.from.bash.to.R[[23]]
here.filePath.output.harmonised.data <- arguments.passed.from.bash.to.R[[24]]
here.filePath.output.MR.analysis <- arguments.passed.from.bash.to.R[[25]]
here.folderPath.output.MR.report <- arguments.passed.from.bash.to.R[[26]]
here.filePath.output.hori.pleio <- arguments.passed.from.bash.to.R[[27]]

# Inspect variables here that are passed from bash arguments
## Total arguments=27
print("============================================================================================")
print("===============Inspect variables here that are passed from bash arguments===================")
print("============================================================================================")

print(paste0("here.exposure.file.path=",here.exposure.file.path))
print(paste0("here.exposure.consortium=",here.exposure.consortium))
print(paste0("here.exposure.trait=",here.exposure.trait))
print(paste0("here.exposure.file.delimiter=",here.exposure.file.delimiter))
print(paste0("here.exposure.file.colname.for.SNP=", here.exposure.file.colname.for.SNP))
print(paste0("here.exposure.file.colname.for.beta=",here.exposure.file.colname.for.beta))
print(paste0("here.exposure.file.colname.for.standard.error=",here.exposure.file.colname.for.standard.error))
print(paste0("here.exposure.file.colname.for.effect.allele=",here.exposure.file.colname.for.effect.allele))
print(paste0("here.exposure.file.colname.for.other.allele=",here.exposure.file.colname.for.other.allele))
print(paste0("here.exposure.file.colname.for.P.value=",here.exposure.file.colname.for.P.value))
print(paste0("here.exposure.clumping.p1=",here.exposure.clumping.p1))
print(paste0("here.exposure.clumping.p2=",here.exposure.clumping.p2))

print(paste0("here.outcome.file.path=",here.outcome.file.path))
print(paste0("here.outcome.consortium=",here.outcome.consortium))
print(paste0("here.outcome.trait=",here.outcome.trait))
print(paste0("here.outcome.file.delimiter=",here.outcome.file.delimiter))
print(paste0("here.outcome.file.colname.for.SNP=",here.outcome.file.colname.for.SNP))
print(paste0("here.outcome.file.colname.for.beta=",here.outcome.file.colname.for.beta))
print(paste0("here.outcome.file.colname.for.standard.error=",here.outcome.file.colname.for.standard.error))
print(paste0("here.outcome.file.colname.for.effect.allele=",here.outcome.file.colname.for.effect.allele))
print(paste0("here.outcome.file.colname.for.other.allele=",here.outcome.file.colname.for.other.allele))
print(paste0("here.outcome.file.colname.for.P.value=",here.outcome.file.colname.for.P.value))

print(paste0("here.iteration=",here.iteration))
print(paste0("here.filePath.output.harmonised.data=",here.filePath.output.harmonised.data))
print(paste0("here.filePath.output.MR.analysis=",here.filePath.output.MR.analysis))
print(paste0("here.folderPath.output.MR.report=",here.folderPath.output.MR.report))
print(paste0("here.filePath.output.hori.pleio=",here.filePath.output.hori.pleio))

# Set up directory under home folder
homeDir <- "/mnt/backedup/home/lunC/"
locRFunction <- paste0(homeDir,"scripts/RFunctions/");

# Set up directory under working
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"

# Import functions
source(paste0(locRFunction,"Rfunction_harmonise-exposures-and-outcome_two-sample-MR.R"))

#----------------------------------------------------------------------------------------------------
# Read a clumped GWAS of trait 1 as an exposure, QCed GWAS of trait 2 as an outcome, harmonise the exposure and outcome, and export the harmonised data set as a tsv file
#----------------------------------------------------------------------------------------------------

# Run analysis
cat("\n"
    ,"===================================================================================================="
    ,"\n", "*************************** Iteration ",here.iteration, "*******************************" 
    ,"\n", "Harmonising exposure ",here.exposure.consortium,"-",here.exposure.trait," and outcome ", here.outcome.consortium,"-", here.outcome.trait 
    ,"\n","==============================================================================================="
    ,"\n")

# Harmonise an exposure and outcome using an external function
## Number of arguments= 193-168+1 = 26
Run.two.sample.Mendelian.randomisation(   exposure.consortium=here.exposure.consortium
                                         ,exposure.trait=here.exposure.trait
                                         ,exposure.file.path=here.exposure.file.path
                                         ,exposure.file.delimiter=here.exposure.file.delimiter
                                         ,exposure.file.colname.for.SNP=here.exposure.file.colname.for.SNP
                                         ,exposure.file.colname.for.beta=here.exposure.file.colname.for.beta
                                         ,exposure.file.colname.for.standard.error=here.exposure.file.colname.for.standard.error
                                         ,exposure.file.colname.for.effect.allele=here.exposure.file.colname.for.effect.allele
                                         ,exposure.file.colname.for.other.allele=here.exposure.file.colname.for.other.allele
                                         ,exposure.file.colname.for.P.value=here.exposure.file.colname.for.P.value
                                         ,exposure.clumping.p1=here.exposure.clumping.p1
                                         ,exposure.clumping.p2=here.exposure.clumping.p2
                                         ,outcome.consortium=here.outcome.consortium
                                         ,outcome.trait=here.outcome.trait
                                         ,outcome.file.path=here.outcome.file.path
                                         ,outcome.file.delimiter=here.outcome.file.delimiter
                                         ,outcome.file.colname.for.SNP=here.outcome.file.colname.for.SNP
                                         ,outcome.file.colname.for.beta=here.outcome.file.colname.for.beta
                                         ,outcome.file.colname.for.standard.error=here.outcome.file.colname.for.standard.error
                                         ,outcome.file.colname.for.effect.allele=here.outcome.file.colname.for.effect.allele
                                         ,outcome.file.colname.for.other.allele=here.outcome.file.colname.for.other.allele
                                         ,outcome.file.colname.for.P.value=here.outcome.file.colname.for.P.value
                                         ,filePath.output.harmonised.data=here.filePath.output.harmonised.data
                                         ,filePath.output.MR.analysis=here.filePath.output.MR.analysis
                                         ,folderPath.output.MR.reports=here.folderPath.output.MR.report
                                         ,filePath.output.horizontal.pleiotropy.test=here.filePath.output.hori.pleio)

#file.copy(paste0(locScripts,"MR_step06-03-02_jobScript-harmonise-exposure-outcome.R"),paste0(locScripts,"MR_step06-03-05_tabulate-two-sample-MR-analysis-results.R"))

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#