#!/usr/bin/env Rscript

##################################################################################
# filename: MR_step08-01-01_create-files-to-loop-through-running-LDSC-munge.R
# Modified from: MR_step06-03-01_create-iterators-for-running-two-sample-MR-on-HPC.R
# program author: Chang
# purpose: 
# date created: 20190410
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Outpu paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20191129  Exported paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv")
# 20191001  Added trait.definitions and data.type. Exported paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv")
# 20190812  Exported the 1 files above
# 20190410  Exported the 1 files above
#-----------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

#-------------------------------------------
# Folder locations under lunC working
#-------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";

locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
locICC <- paste0(locMR,"ICC-cannabis-ever/")

locGSCAN.QC3 <- paste0(locMR,"noICC_results/QC3_remove_ambiguousSNPs_indel")

locUKB.3456.QC3 <- paste0(locMR,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")

locUKB.ESDPW.QC3 <- paste0(locMR,"UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")

locUKB.CCPD.QC3 <- paste0(locMR,"UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")

locUKB.20161.QC3 <- paste0(locMR,"UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")

locUKB.caffeine.QC3 <- paste0(locMR,"UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")

# Output folders
loc.LDSC <- paste0(workingDir,"MR_ICC_GSCAN_201806/LD-score-correlation/")
loc.LDSC.input <- paste0(loc.LDSC,"input/")

dir.create(loc.LDSC)
dir.create(loc.LDSC.input)

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#------------------------------------------------------------------------------
# Get file paths of QCed GWASs
#------------------------------------------------------------------------------
filePath.QCed.GWAS.GSCAN <- list.files(path=c(locGSCAN.QC3)
                               ,pattern = "*_noICC.ambiguousSNPRemoved$"
                               ,full.names = TRUE) # length(filePath.QCed.GWAS.GSCAN) 5

filePath.QCed.GWAS.ICC <- Sys.glob(paste0(locICC,"Cannabis_ICC_UKB_small.txt")) # length(filePath.QCed.GWAS.ICC) 1

filePath.QCed.GWAS.UKB <- list.files(path=c( locUKB.3456.QC3
                                            ,locUKB.ESDPW.QC3
                                            ,locUKB.CCPD.QC3
                                            ,locUKB.20161.QC3
                                            ,locUKB.caffeine.QC3)
                                     ,pattern= glob2rx("^QCed-GWAS-UKB*_headed$")
                                     ,full.names = TRUE) # length(filePath.QCed.GWAS.UKB) 5

#----------------------------------------------------------------------------------------------------
# Add information needed for LDSC to a file
#----------------------------------------------------------------------------------------------------
# Directory of QCed GWAS files
file.path.QCed.GWASs <- c( filePath.QCed.GWAS.GSCAN
                           ,filePath.QCed.GWAS.ICC
                           ,filePath.QCed.GWAS.UKB) # length(file.path.QCed.GWASs) 11

# Create outcome file information
info.QCed.GWASs <- data.frame(filePath=file.path.QCed.GWASs
                              ,colname.SNP=c(rep("RSID",times=length(filePath.QCed.GWAS.GSCAN))
                                             ,rep("SNP",times=length(filePath.QCed.GWAS.ICC))
                                             ,rep("SNP",times=length(filePath.QCed.GWAS.UKB)))
                              ,colname.effect.allele=c(rep("ALT",times=length(filePath.QCed.GWAS.GSCAN))
                                                       ,rep("Allele1",times=length(filePath.QCed.GWAS.ICC))
                                                       ,rep("ALLELE1",times=length(filePath.QCed.GWAS.UKB)))
                              ,colname.other.allele=c(rep("REF",times=length(filePath.QCed.GWAS.GSCAN))
                                                      ,rep("Allele2",times=length(filePath.QCed.GWAS.ICC))
                                                      ,rep("ALLELE0",times=length(filePath.QCed.GWAS.UKB)))
                              ,colname.beta=c(rep("BETA",times=length(filePath.QCed.GWAS.GSCAN))
                                              ,rep("Effect",times=length(filePath.QCed.GWAS.ICC))
                                              ,rep("BETA",times=length(filePath.QCed.GWAS.UKB)))
                              ,colname.SE=c(rep("SE",times=length(filePath.QCed.GWAS.GSCAN))
                                            ,rep("StdErr",times=length(filePath.QCed.GWAS.ICC))
                                            ,rep("SE",times=length(filePath.QCed.GWAS.UKB)))
                              ,colname.p.value=c(rep("PVALUE",times=length(filePath.QCed.GWAS.GSCAN))
                                                 ,rep("P",times=length(filePath.QCed.GWAS.ICC))
                                                 ,rep("PVALUE",times=length(filePath.QCed.GWAS.UKB)))
                              ,consortium=c(rep("GSCAN",times=length(filePath.QCed.GWAS.GSCAN))
                                            ,rep("ICC",times=length(filePath.QCed.GWAS.ICC))
                                            ,rep("UKB",times=length(filePath.QCed.GWAS.UKB)))
                              ,trait=c( "AI","CPD","DPW","SC","SI"
                                        ,"CI"
                                        ,"CCPD","caffeine","ESDPW","PYOS","CPD")
                              ,trait.type=c("continuous","continuous","continuous","binary","binary"
                                            ,"binary"
                                            ,"continuous","continuous","continuous","continuous","continuous")
                              ,N.continuous=c(119239,122027,185828,NA,NA
                                              ,NA
                                              ,276533,168919,296735,186411,25153)
                              ,N.cases=c(NA,NA,NA,82738,112172
                                         ,44084
                                         ,NA,NA,NA,NA,NA)
                              ,N.controls=c(NA,NA,NA,42623,95553.96
                                            ,120657
                                            ,NA,NA,NA,NA,NA)
                              ,population.prevalence=c(NA,NA,NA,0.25,0.45
                                                       ,0.035
                                                       ,NA,NA,NA,NA,NA)
                              ,file.delimiter=c(rep("tab",times=length(filePath.QCed.GWAS.GSCAN))
                                                ,rep("space",times=length(filePath.QCed.GWAS.ICC))
                                                ,rep("tab",times=length(filePath.QCed.GWAS.UKB)))
                              # Add trait definitions
                              ,trait.definition=c("Age of initiating regular smoking","Cigarettes smoked per day using averaged number or 5 response categories: 1-5, 6-15, 16-25, 26-35, or 36+","Number of drinks consumed per week in ever drinkers","Current smokers versus former smokers","Ever versus never smoked regularly. A regular smoker was defined as having smoked > 100 cigarettes during lifetime, ever having smoked every day for at least one month, or simply ever smoking regularly","Ever versus never taken cannabis","Number of cups of coffee consumed per day (include decaffeinated coffee; UKB data field: 1498)","Caffeine consumed per day through regular coffee and tea (mg/day; UKB Data-Field: 1488, 1498, 1508)","Number of standard drinks consumed per week based on alcohol intake frequency, types of alcoholic beverage and intake quantity","Pack years of smoking in ever or current cigarette smokers (UKB Data-Field: 20161)","Number of cigarettes currently smoked daily (UKB Data-Field: 3456)")
                              ,data.type=c(rep("GWAS_meta-analysis",times=5)
                                           ,"GWAS"
                                           ,rep("GWAS",times=5))
                              ,substance=dplyr::case_when( trait %in% c("AI","CPD","PYOS","SC","SI") ~ "nicotine"
                                                          ,trait %in% c("DPW","ESDPW") ~ "alcohol"
                                                          ,trait %in% c("CI") ~ "cannabis"
                                                          ,trait %in% c("caffeine") ~ "caffeine"
                                                          ,TRUE ~ as.character(trait))
                              ,stringsAsFactors = FALSE) # dim(info.QCed.GWASs) 11 18

# Export files
ExportFileTabSeparated( data=info.QCed.GWASs
                        ,missing.values.as = "nan"
                        ,output.file.path = paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv"))

#file.copy(paste0(locScripts,"MR_step08-01-01_create-files-to-loop-through-running-LDSC-munge.R"),paste0(locScripts,"MR_step08-02-01_create-files-to-loop-through-running-LDSC-genetic-correlations.R"))

#-------------------------------------------------------------------------------------------------
#------------------This is the end of this file---------------------------------------------------
#-------------------------------------------------------------------------------------------------
