#####################################################################################################################
# filename: MR_step06-02_prepare-input-files-for-twoSampleMR-analysis.R
# program author: Chang
# purpose: Process exposure GWAS files for two sample MR analysis : (1)add sample sizes to UKB clumped GWAS, (2) standardise file format (all space-separated) for running two-sample-MR. GSCAN and ICC GWASs came with sample size columns.
# date created: 20190307
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Input paste0(locICCRef,"top-hits-substance-use.txt")

# Input	paste0(locGSCAN,"GWAS_from-clumped-SNPs_dpw_noICC_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
# Input	paste0(locGSCAN,"GWAS_from-clumped-SNPs_dpw_noICC_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")

# Input	paste0(ref.UKB.3456CPD,"BOLT-LMM-ukb3456_IID_NA_in_20453-pheno-X3456_mean/bolt_lmm_ukb_imp_chr1_v3.bgen-X3456_mean.log")
# Input	paste0(locUKB.3456.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
# Input	paste0(locUKB.3456.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")

# Input	paste0(ref_UKB_ESDPW,"BOLT-LMM-phenotype-pheno-complete_alcohol_unitsweekly/bolt_lmm_ukb_imp_chr9_v3.bgen-complete_alcohol_unitsweekly.log")
# Input	paste0(locUKB.ESDPW.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
# Input	paste0(locUKB.ESDPW.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")
# Input	paste0(ref_UKB_CCPD,"BOLT-LMM-phenotype-pheno-all_coffee_cpd/bolt_lmm_ukb_imp_chr22_v3.bgen-all_coffee_cpd.log")
# Input	paste0(locUKB.CCPD.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
# Input	paste0(locUKB.CCPD.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")

# Input	paste0(ref_UKB20161,"BOLT-LMM-phenotype-pheno-merged_pack_years_20161/bolt_lmm_ukb_imp_chr9_v3.bgen-merged_pack_years_20161.log")
# Input	paste0(locUKB.20161.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
# Input	paste0(locUKB.20161.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")

# Outpu paste0(locGSCAN,nam.clumped.GWAS.GSCAN.DPW.p1.5e8,"_linear-BETA-added")
# Outpu paste0(locGSCAN,nam.clumped.GWAS.GSCAN.DPW.p1.1e5,"_linear-BETA-added")

# Outpu paste0(locUKB.3456.QC4,nam.clumped.GWAS.UKB.3456CPD.p1.5e8,"_sample-size-added")
# Outpu paste0(locUKB.3456.QC4,nam.clumped.GWAS.UKB.3456CPD.p1.1e5,"_sample-size-added")

# Outpu paste0(locUKB.ESDPW.QC4,nam.clumped.GWAS.UKB.ESDPW.p1.5e8,"_sample-size-added")
# Outpu paste0(locUKB.ESDPW.QC4,nam.clumped.GWAS.UKB.ESDPW.p1.1e5,"_sample-size-added")

# Outpu paste0(locUKB.CCPD.QC4,nam.clumped.GWAS.UKB.CCPD.p1.5e8,"_sample-size-added")
# Outpu paste0(locUKB.CCPD.QC4,nam.clumped.GWAS.UKB.CCPD.p1.1e5,"_sample-size-added")

# Outpu paste0(locUKB.caffeine.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5","_sample-size-added")
# Outpu paste0(locUKB.caffeine.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6","_sample-size-added")

# Outpu paste0(locUKB.20161.QC4,nam.clumped.GWAS.UKB.PYOS.p1.5e8,"_sample-size-added")
# Outpu paste0(locUKB.20161.QC4,nam.clumped.GWAS.UKB.PYOS.p1.1e5,"_sample-size-added")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190812  Exported 2 files for UKB caffeine (all space-separated)
# 20190405  Exported the 10 files above (all space-separated)
# 20190307  Exported 9 files above (all space-separated)
#-----------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locPlots_MR <- paste0(homeDir,"plots/MR_ICC_GSCAN_201806/")

#-------------------------------------------
# Folder locations under Stuartma lab
#-------------------------------------------
stuartmaDir <- "/reference/data/UKBB_500k/versions/lab_stuartma/"
stuartma_pheno <- paste0(stuartmaDir,"pheno/")

ref_UKB_GWAS_BOLT_LMM <- paste0(stuartmaDir,"draft_gwas/BOLT_LMM/")
ref_UKB_ESDPW <- paste0(ref_UKB_GWAS_BOLT_LMM,"UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/")
ref_UKB_CCPD <- paste0(ref_UKB_GWAS_BOLT_LMM,"UKB-cups-of-coffee-per-day_IID-NA-in-UKB20453-everUsedCannabis/")
ref_UKB20161 <- paste0(ref_UKB_GWAS_BOLT_LMM,"UKB20161-pack-years-of-smoking_IID-NA-in-UKB20453-everUsedCannabis/")
ref.UKB.3456CPD <- paste0(ref_UKB_GWAS_BOLT_LMM,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/")
ref.UKB.caffeine <- paste0(ref_UKB_GWAS_BOLT_LMM,"UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/")

# Folders under lunC home  
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
locUKB3456_pheno <- paste0(homeDir,"data/UKBionbank_phenotype/ukb3456_numCigareDaily/")

# Folders under lunC working
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locICCRef <- paste0(workingDir,"MR_ICC_GSCAN_201806/reference/")

locICC <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
locGSCAN <- paste0(locICC,"noICC_results/QC4_GWAS_from_clumped_SNPs/")

#--------------------------------------------------------------------
# Folder locations under my working 
#--------------------------------------------------------------------
# UKB GWAS 
locUKB.3456.QC3 <- paste0(locICC,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/")
locUKB.3456.QC4 <- paste0(locICC,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")

locUKB.ESDPW.QC3 <- paste0(locICC,"UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/")
locUKB.ESDPW.QC4 <- paste0(locICC,"UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")


locUKB.CCPD.QC3 <- paste0(locICC,"UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/")
locUKB.CCPD.QC4 <- paste0(locICC,"UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")

locUKB.20161.QC3 <- paste0(locICC,"UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/")
locUKB.20161.QC4 <- paste0(locICC,"UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")

locUKB.caffeine.QC3 <- paste0(locICC,"UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/")
locUKB.caffeine.QC4 <- paste0(locICC,"UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")

locMR_multiSNPs <- paste0(workingDir,"MR_ICC_GSCAN_201806/result_MR_multipleSNPs/")

library(dplyr, lib.loc = "/software/R/R-3.4.1/lib64/R/library")
library(stringr, lib.loc = "/software/R/R-3.4.1/lib64/R/library")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#--------------------------------------------------------------------------------------------------
# --------------------Add linear beta to GSCAN clumped GWAS of drinks per week
#--------------------------------------------------------------------------------------------------
nam.clumped.GWAS.GSCAN.DPW.p1.5e8 <- "GWAS_from-clumped-SNPs_dpw_noICC_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
nam.clumped.GWAS.GSCAN.DPW.p1.1e5 <- "GWAS_from-clumped-SNPs_dpw_noICC_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5"
  
ImportASpaceSeparatedFile(input.file.path = paste0(locGSCAN,nam.clumped.GWAS.GSCAN.DPW.p1.5e8)
                            , data.name = "dat.clumped.GWAS.GSCAN.DPW.p1.5e8") # dim(dat.clumped.GWAS.GSCAN.DPW.p1.5e8) 5 14

ImportASpaceSeparatedFile(input.file.path = paste0(locGSCAN,nam.clumped.GWAS.GSCAN.DPW.p1.1e5)
                          , data.name = "dat.clumped.GWAS.GSCAN.DPW.p1.1e5") # dim(dat.clumped.GWAS.GSCAN.DPW.p1.1e5) 82 14

# Compute linear betas for each DPW clumped SNP. Because the phenotype drinks per week was log-transformed, the BETAs in GSCAN GWAS files are in log scale, which is difficul to interpret.   

## step 1: Copy "Mean(Var)" under Drinks per Week from Supplementary Table 1
GSCAN_DPW_meanVar <- c(1.3,1.2,1.2,NA,0.8,2.5,NA,1.6,0.6,0.2,0.3,1.6,2.0,7.3,1.2,1.2,0.9,1.1,2.9,0.5,0.5,2.0,1.6,2.1,1.8 ,1.9,1.9,2.2,2.3,2.2,0.6,0.6,0.6,0.7,0.7,0.7)

## step 2: Take an average
GSCAN_DPW_meanVar_avg <- mean(GSCAN_DPW_meanVar,na.rm = T)

## step 3: For each clumped SNP, calculate its linear beta using this equation: exp(mean in log scale+ beta in log scale)- exp(mean in log scale)= linear beta
dat.clumped.GWAS.GSCAN.DPW.p1.5e8$BETA.linear <- with(dat.clumped.GWAS.GSCAN.DPW.p1.5e8, exp(GSCAN_DPW_meanVar_avg+ BETA)-exp(GSCAN_DPW_meanVar_avg))

dat.clumped.GWAS.GSCAN.DPW.p1.1e5$BETA.linear <- with(dat.clumped.GWAS.GSCAN.DPW.p1.1e5, exp(GSCAN_DPW_meanVar_avg+ BETA)-exp(GSCAN_DPW_meanVar_avg))

# Export data 
ExportFileSpaceSeparated(data = dat.clumped.GWAS.GSCAN.DPW.p1.5e8
                         ,output.file.path = paste0(locGSCAN,nam.clumped.GWAS.GSCAN.DPW.p1.5e8,"_linear-BETA-added"))

ExportFileSpaceSeparated(data = dat.clumped.GWAS.GSCAN.DPW.p1.1e5
                         ,output.file.path = paste0(locGSCAN,nam.clumped.GWAS.GSCAN.DPW.p1.1e5,"_linear-BETA-added"))

#---------------------------------------------------------------------------------------
# Process exposure GWAS file of UKB 3456 excluding ever using cannabis 20453
## SNPs selected via LD clumping
#---------------------------------------------------------------------------------------
# Get sample size from any BOLT-LMM log file 25153

## Import a log file that has no delimiters
log.file <- readLines(paste0(ref.UKB.3456CPD
                            ,"BOLT-LMM-ukb3456_IID_NA_in_20453-pheno-X3456_mean/"
                            ,"bolt_lmm_ukb_imp_chr1_v3.bgen-X3456_mean.log"))
## Extract the line with a pattern specified by regular expression
### double backslash to escape the special character ( and )
pattern.sample.size <- "^Number of indivs with no missing phenotype\\(s\\) to use"
line.extracted <- grep(pattern.sample.size
                       ,log.file
                       ,value = T)

N.overall.UKB3456 <- stringr::str_replace(string=line.extracted
                                          ,pattern = "^Number of indivs with no missing phenotype\\(s\\) to use\\: "
                                          ,replacement = "") # "25153"

## Impport clumped GWAS file
nam.clumped.GWAS.UKB.3456CPD.p1.5e8 <- "GWAS_from-clumped-SNPs_GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
nam.clumped.GWAS.UKB.3456CPD.p1.1e5 <- "GWAS_from-clumped-SNPs_GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5"

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.3456.QC4,nam.clumped.GWAS.UKB.3456CPD.p1.5e8)
                          ,data.name = "dat.clumped.GWAS.UKB.3456CPD.p1.5e8") # dim(dat.clumped.GWAS.UKB.3456CPD.p1.5e8) 3 18

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.3456.QC4,nam.clumped.GWAS.UKB.3456CPD.p1.1e5)
                          ,data.name = "dat.clumped.GWAS.UKB.3456CPD.p1.1e5") # dim(dat.clumped.GWAS.UKB.3456CPD.p1.1e5) 48 18

dat.clumped.GWAS.UKB.3456CPD.p1.5e8$N.overall <- as.numeric(N.overall.UKB3456)
dat.clumped.GWAS.UKB.3456CPD.p1.1e5$N.overall <- as.numeric(N.overall.UKB3456)

# Export the processed file to input folder
ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.3456CPD.p1.5e8
                         ,output.file.path = paste0(locUKB.3456.QC4
                                                    ,nam.clumped.GWAS.UKB.3456CPD.p1.5e8
                                                    ,"_sample-size-added"))

ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.3456CPD.p1.1e5
                         ,output.file.path = paste0(locUKB.3456.QC4
                                                    ,nam.clumped.GWAS.UKB.3456CPD.p1.1e5
                                                    ,"_sample-size-added"))

#---------------------------------------------------------------------------------------
# Process exposure GWAS file of UKB ESDPW excluding ever using cannabis 20453
## SNPs selected via LD clumping
#---------------------------------------------------------------------------------------
# Get sample size from any BOLT-LMM log file 

## Import a log file that has no delimiters
log.file.UKB.ESDPW <- readLines(paste0(ref_UKB_ESDPW
                                       ,"BOLT-LMM-phenotype-pheno-complete_alcohol_unitsweekly/"
                                       ,"bolt_lmm_ukb_imp_chr9_v3.bgen-complete_alcohol_unitsweekly.log"))

## Extract the line with a pattern specified by regular expression
### double backslash to escape the special character ( and )
pattern.sample.size <- "^Number of indivs with no missing phenotype\\(s\\) to use"
log.file.UKB.ESDPW.extracted.line <- grep(pattern.sample.size
                                          ,log.file.UKB.ESDPW
                                          ,value = T)

N.overall.UKB.ESDPW <- stringr::str_replace(string=log.file.UKB.ESDPW.extracted.line
                                          ,pattern = "^Number of indivs with no missing phenotype\\(s\\) to use\\: "
                                          ,replacement = "") # "296735"

## Impport clumped GWAS file
nam.clumped.GWAS.UKB.ESDPW.p1.5e8 <- "GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
nam.clumped.GWAS.UKB.ESDPW.p1.1e5 <- "GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5"
  
ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.ESDPW.QC4,nam.clumped.GWAS.UKB.ESDPW.p1.5e8)
                          ,data.name = "dat.clumped.GWAS.UKB.ESDPW.p1.5e8") # dim(dat.clumped.GWAS.UKB.ESDPW.p1.5e8) 39 18

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.ESDPW.QC4,nam.clumped.GWAS.UKB.ESDPW.p1.1e5)
                          ,data.name = "dat.clumped.GWAS.UKB.ESDPW.p1.1e5") # dim(dat.clumped.GWAS.UKB.ESDPW.p1.1e5) 221  18

dat.clumped.GWAS.UKB.ESDPW.p1.5e8$N.overall <- as.numeric(N.overall.UKB.ESDPW)
dat.clumped.GWAS.UKB.ESDPW.p1.1e5$N.overall <- as.numeric(N.overall.UKB.ESDPW)

# Export the processed file to input folder
ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.ESDPW.p1.5e8  
                         ,output.file.path = paste0(locUKB.ESDPW.QC4
                                                    ,nam.clumped.GWAS.UKB.ESDPW.p1.5e8
                                                    ,"_sample-size-added"))

ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.ESDPW.p1.1e5  
                         ,output.file.path = paste0(locUKB.ESDPW.QC4
                                                    ,nam.clumped.GWAS.UKB.ESDPW.p1.1e5
                                                    ,"_sample-size-added"))

#---------------------------------------------------------------------------------------
# Process exposure GWAS file of UKB CCPD excluding ever using cannabis 20453
## SNPs selected via LD clumping
#---------------------------------------------------------------------------------------
# Get sample size from any BOLT-LMM log file 

## Import a log file that has no delimiters
log.file.UKB.CCPD <- readLines(paste0(ref_UKB_CCPD
                                       ,"BOLT-LMM-phenotype-pheno-all_coffee_cpd/"
                                       ,"bolt_lmm_ukb_imp_chr22_v3.bgen-all_coffee_cpd.log"))

## Extract the line with a pattern specified by regular expression
### double backslash to escape the special character ( and )
pattern.sample.size <- "^Number of indivs with no missing phenotype\\(s\\) to use"
log.file.UKB.CCPD.extracted.line <- grep(pattern.sample.size
                                          ,log.file.UKB.CCPD
                                          ,value = T)

N.overall.UKB.CCPD <- stringr::str_replace(string=log.file.UKB.CCPD.extracted.line
                                            ,pattern = "^Number of indivs with no missing phenotype\\(s\\) to use\\: "
                                            ,replacement = "") # "276533"

## Impport clumped GWAS file
nam.clumped.GWAS.UKB.CCPD.p1.5e8 <- "GWAS_from-clumped-SNPs_GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
nam.clumped.GWAS.UKB.CCPD.p1.1e5 <- "GWAS_from-clumped-SNPs_GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5"

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.CCPD.QC4,nam.clumped.GWAS.UKB.CCPD.p1.5e8)
                          ,data.name = "dat.clumped.GWAS.UKB.CCPD.p1.5e8") # dim(dat.clumped.GWAS.UKB.CCPD.p1.5e8) 25 17

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.CCPD.QC4,nam.clumped.GWAS.UKB.CCPD.p1.1e5)
                          ,data.name = "dat.clumped.GWAS.UKB.CCPD.p1.1e5") # dim(dat.clumped.GWAS.UKB.CCPD.p1.1e5) 110  17

dat.clumped.GWAS.UKB.CCPD.p1.5e8$N.overall <- as.numeric(N.overall.UKB.CCPD)
dat.clumped.GWAS.UKB.CCPD.p1.1e5$N.overall <- as.numeric(N.overall.UKB.CCPD)

# Export the processed file to input folder
ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.CCPD.p1.5e8
                         ,output.file.path = paste0(locUKB.CCPD.QC4
                                                    ,nam.clumped.GWAS.UKB.CCPD.p1.5e8
                                                    ,"_sample-size-added"))

ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.CCPD.p1.1e5
                         ,output.file.path = paste0(locUKB.CCPD.QC4
                                                    ,nam.clumped.GWAS.UKB.CCPD.p1.1e5
                                                    ,"_sample-size-added"))

#---------------------------------------------------------------------------------------
# Process exposure GWAS file of UKB caffeine consumed per day excluding ever using cannabis 20453
## SNPs selected via LD clumping
#---------------------------------------------------------------------------------------
# Get sample size from any BOLT-LMM log file 

## Import a log file that has no delimiters
## Extract the line with a pattern specified by regular expression
### double backslash to escape the special character ( and )

pattern.sample.size <- "^Number of indivs with no missing phenotype\\(s\\) to use\\: "
N.overall.UKB.caffeine <- readLines(paste0(ref.UKB.caffeine
                                           ,"BOLT-LMM-phenotype-pheno-caffeine.per.day/"
                                           ,"bolt_lmm_ukb_imp_chr1_v3.bgen-caffeine.per.day.log")) %>% 
  grep(pattern=pattern.sample.size,value=T) %>%
  gsub(pattern=pattern.sample.size, replacement="")

# Import clumped GWAS files and add sample size
dat.clumped.GWAS.UKB.caffeine.p1.1e5 <- ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.caffeine.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")
                                                                  ,data.name = "dat.clumped.GWAS.UKB.caffeine.p1.1e5") %>%
  dplyr::mutate(N.overall=as.numeric(N.overall.UKB.caffeine)) # dim(dat.clumped.GWAS.UKB.caffeine.p1.1e5) 103 18

dat.clumped.GWAS.UKB.caffeine.p1.5e8 <- ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.caffeine.QC4,"GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
                                                                  ,data.name = "dat.clumped.GWAS.UKB.caffeine.p1.5e8") %>%
  dplyr::mutate(N.overall=as.numeric(N.overall.UKB.caffeine)) # dim(dat.clumped.GWAS.UKB.caffeine.p1.5e8) 22 18

# Export the processed file to input folder
ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.caffeine.p1.1e5
                         ,output.file.path = paste0(locUKB.caffeine.QC4
                                                    ,"GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5"
                                                    ,"_sample-size-added"))

ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.caffeine.p1.5e8
                         ,output.file.path = paste0(locUKB.caffeine.QC4
                                                    ,"GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
                                                    ,"_sample-size-added"))

#---------------------------------------------------------------------------------------
# Process exposure GWAS file of UKB PYOS  excluding ever using cannabis 20453
## SNPs selected via LD clumping
#---------------------------------------------------------------------------------------
# Get sample size from any BOLT-LMM log file 

## Import a log file that has no delimiters
log.file.UKB.PYOS <- readLines(paste0(ref_UKB20161
                                      ,"BOLT-LMM-phenotype-pheno-merged_pack_years_20161/"
                                      ,"bolt_lmm_ukb_imp_chr9_v3.bgen-merged_pack_years_20161.log"))

## Extract the line with a pattern specified by regular expression
### double backslash to escape the special character ( and )
pattern.sample.size <- "^Number of indivs with no missing phenotype\\(s\\) to use\\: "
log.file.UKB.PYOS.extracted.line <- grep(pattern.sample.size
                                         ,log.file.UKB.PYOS
                                         ,value = T)

N.overall.UKB.PYOS <- stringr::str_replace(string=log.file.UKB.PYOS.extracted.line
                                           ,pattern = pattern.sample.size
                                           ,replacement = "") # "186411"

## Impport clumped GWAS file
nam.clumped.GWAS.UKB.PYOS.p1.5e8 <- "GWAS_from-clumped-SNPs_GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
nam.clumped.GWAS.UKB.PYOS.p1.1e5 <- "GWAS_from-clumped-SNPs_GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5"

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.20161.QC4,nam.clumped.GWAS.UKB.PYOS.p1.5e8)
                          ,data.name = "dat.clumped.GWAS.UKB.PYOS.p1.5e8") # dim(dat.clumped.GWAS.UKB.PYOS.p1.5e8) 20 17

ImportASpaceSeparatedFile(input.file.path = paste0(locUKB.20161.QC4,nam.clumped.GWAS.UKB.PYOS.p1.1e5)
                          ,data.name = "dat.clumped.GWAS.UKB.PYOS.p1.1e5") # dim(dat.clumped.GWAS.UKB.PYOS.p1.1e5) 160  17

dat.clumped.GWAS.UKB.PYOS.p1.5e8$N.overall <- as.numeric(N.overall.UKB.PYOS)
dat.clumped.GWAS.UKB.PYOS.p1.1e5$N.overall <- as.numeric(N.overall.UKB.PYOS)

# Export the processed file to input folder
ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.PYOS.p1.5e8
                         ,output.file.path = paste0(locUKB.20161.QC4
                                                    ,nam.clumped.GWAS.UKB.PYOS.p1.5e8
                                                    ,"_sample-size-added"))

ExportFileSpaceSeparated(data= dat.clumped.GWAS.UKB.PYOS.p1.1e5
                         ,output.file.path = paste0(locUKB.20161.QC4
                                                    ,nam.clumped.GWAS.UKB.PYOS.p1.1e5
                                                    ,"_sample-size-added"))

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#
