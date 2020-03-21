#!/usr/bin/env Rscript

##################################################################################
# File name: MR_step06-03-01_create-iterators-for-running-two-sample-MR-on-HPC.R
# Old file name: zMR_step06-03_two-sample-MR.R
# program author: Chang
# purpose: Analyse causal relationship between exposure GWAS (QCed and clumped/COJOed) and outcome GWAS (QCed)
# date created: 20190405
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Outpu paste0(loc.twoSampleMR.input,"filePath_exposure-clumped-GWASs.tsv")
# Outpu paste0(loc.twoSampleMR.input,"filePath_outcome-QCed-GWASs.tsv")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190812  Exported the 2 files above
# 20190408  Exported the 2 files above
#-----------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

#-------------------------------------------
# Folders under lunC working directory
#-------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data

locICC <- paste0(locMR,"ICC-cannabis-ever/")
locICC.QC4 <- paste0(locMR,"ICC-cannabis-ever/QC4_GWAS_from_clumped_SNPs")

locGSCAN.QC3 <- paste0(locMR,"noICC_results/QC3_remove_ambiguousSNPs_indel")
locGSCAN.QC4 <- paste0(locMR,"noICC_results/QC4_GWAS_from_clumped_SNPs") # location of exposure data files

locUKB.3456.QC3 <- paste0(locMR,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")
locUKB.3456.QC4 <- paste0(locMR,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs")

locUKB.ESDPW.QC3 <- paste0(locMR,"UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")
locUKB.ESDPW.QC4 <- paste0(locMR,"UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs")

locUKB.CCPD.QC3 <- paste0(locMR,"UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")
locUKB.CCPD.QC4 <- paste0(locMR,"UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs")

locUKB.20161.QC3 <- paste0(locMR,"UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")
locUKB.20161.QC4 <- paste0(locMR,"UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs")

locUKB.caffeine.QC3 <- paste0(locMR,"UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC3_remove_ambiguousSNPs_indel")
locUKB.caffeine.QC4 <- paste0(locMR,"UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC4_GWAS_from_clumped_SNPs")

#--------------------------------------------------------------------
# Folder locations for output files
#--------------------------------------------------------------------
loc.twoSampleMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/two-sample-MR/")
loc.twoSampleMR.input <- paste0(loc.twoSampleMR,"input/")
#dir.create(loc.twoSampleMR)
#dir.create(loc.twoSampleMR.input)

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#---------------------------------------------------------------------
# Get file paths of UKB clumped GWAS to use as exposures
#---------------------------------------------------------------------
filePath.clumped.GWAS.UKB <- list.files(path=c(locUKB.3456.QC4
                                               ,locUKB.ESDPW.QC4
                                               ,locUKB.CCPD.QC4
                                               ,locUKB.20161.QC4
                                               ,locUKB.caffeine.QC4)
                                        ,pattern = glob2rx("^GWAS_from-clumped-SNPs_GWAS-UKB*_sample-size-added$")
                                        ,full.names = T) # length(filePath.clumped.GWAS.UKB) 10

#---------------------------------------------------------------------
# Get file paths of GSCAN clumped GWAS to use as exposures
#---------------------------------------------------------------------
pattern.GSCAN.files <- "^GWAS_from-clumped-SNPs_ai_noICC_LDWindow-kb-10000_R2-0.01_*|^GWAS_from-clumped-SNPs_cpd_noICC_LDWindow-kb-10000_R2-0.01*|^GWAS_from-clumped-SNPs_dpw_*linear-BETA-added$|^GWAS_from-clumped-SNPs_sc_noICC_LDWindow-kb-10000_R2-0.01_*|^GWAS_from-clumped-SNPs_si_noICC_LDWindow-kb-10000_R2-0.01_*"

filePath.clumped.GWAS.GSCAN <- list.files(path=locGSCAN.QC4
                                          ,pattern=glob2rx(pattern.GSCAN.files)
                                          ,full.names = TRUE) # length(filePath.clumped.GWAS.GSCAN) 10

#------------------------------------------------------------------------------
# Get file paths of ICC clumped GWAS of cannabis initiation to use as exposures
#------------------------------------------------------------------------------
filePath.clumped.GWAS.ICC <- list.files(path=locICC.QC4
                                        ,pattern="^GWAS_from-clumped-SNPs_GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01"
                                        ,full.names = TRUE) # length(filePath.clumped.GWAS.ICC) 2

#------------------------------------------------------------------------------
# Get file paths of QCed GWASs to use as outcomes
#------------------------------------------------------------------------------
filePath.QCed.GWAS.GSCAN <- list.files(path=c(locGSCAN.QC3)
                               ,pattern = "*_noICC.ambiguousSNPRemoved$"
                               ,full.names = TRUE) # length(filePath.QCed.GWAS.GSCAN) 5

filePath.QCed.GWAS.ICC <- Sys.glob(paste0(locICC,"Cannabis_ICC_UKB_small.txt")) # length(filePath.QCed.GWAS.ICC) 1

# Note that output will be sorted 
filePath.QCed.GWAS.UKB <- list.files(path=c(locUKB.3456.QC3
                                            ,locUKB.ESDPW.QC3
                                            ,locUKB.CCPD.QC3
                                            ,locUKB.20161.QC3
                                            ,locUKB.caffeine.QC3)
                                     ,pattern= glob2rx("^QCed-GWAS-UKB*_headed$")
                                     ,full.names = TRUE) # length(filePath.QCed.GWAS.UKB) 5

#----------------------------------------------------------------------------------------------------
# Create exposure file and outcome file information for running two sample MR
## MR is run on every exposure-outcome pair, except for exposure-outcome of same consortium and trait
#----------------------------------------------------------------------------------------------------
# Grep clumping criteria from exposure file names
# Search variables to reshape
patterns_to_search <- glob2rx("LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6|LDWindow-kb-10000_R2-0.01_p1-1e-5_p2-1e-5")
exposure.file.paths <- c( filePath.clumped.GWAS.GSCAN
                          ,filePath.clumped.GWAS.ICC
                          ,filePath.clumped.GWAS.UKB) # length(exposure.file.paths) 22

tem1 <- gsub(exposure.file.paths,pattern = "_sample-size-added",replacement = "")
tem2 <- gsub(tem1,pattern="GWAS_from-clumped-SNPs_",replacement="")
tem3 <- gsub(basename(tem2),pattern = "_", replacement = "-")

exposure.file.name.selective <- paste0("clumped-",tem3) # length(exposure.file.name.selective) 22

exposure.GWAS <- data.frame(filePath= exposure.file.paths 
                            ,fileNameSpecial=exposure.file.name.selective
                            ,consortium=c(rep("GSCAN",times=length(filePath.clumped.GWAS.GSCAN))
                                          ,rep("ICC",times=length(filePath.clumped.GWAS.ICC))
                                          ,rep("UKB",times=length(filePath.clumped.GWAS.UKB)))
                            ,trait=c( rep(c("AI","CPD","DPW","SC","SI"),each=2)
                                     ,rep("CI",times=2)
                                     ,rep(c("CCPD","caffeine","ESDPW","PYOS","CPD"),each=2))
                            ,delimiter=rep("space",times=length(exposure.file.paths))
                            ,colname.for.SNP=c(rep("SNP",times=length(exposure.file.paths)))
                            ,colname.for.beta=c(rep("BETA",times=length(filePath.clumped.GWAS.GSCAN))
                                                ,rep("Effect",times=length(filePath.clumped.GWAS.ICC))
                                                ,rep("BETA",times=length(filePath.clumped.GWAS.UKB)))
                            ,colname.for.SE=c(rep("SE",times=length(filePath.clumped.GWAS.GSCAN))
                                              ,rep("StdErr",times=length(filePath.clumped.GWAS.ICC))
                                              ,rep("SE",times=length(filePath.clumped.GWAS.UKB))) 
                            ,colname.for.effect.allele=c(rep("ALT",times=length(filePath.clumped.GWAS.GSCAN))
                                                         ,rep("Allele1",times=length(filePath.clumped.GWAS.ICC))
                                                         ,rep("ALLELE1",times=length(filePath.clumped.GWAS.UKB)))  
                            ,colname.for.other.allele=c(rep("REF",times=length(filePath.clumped.GWAS.GSCAN))
                                                        ,rep("Allele2",times=length(filePath.clumped.GWAS.ICC))
                                                        ,rep("ALLELE0",times=length(filePath.clumped.GWAS.UKB)))
                            ,colname.for.P.value= c(rep("PVALUE",times=length(filePath.clumped.GWAS.GSCAN))
                                                        ,rep("P",times=length(filePath.clumped.GWAS.ICC))
                                                        ,rep("PVALUE",times=length(filePath.clumped.GWAS.UKB)))
                            ,stringsAsFactors = FALSE) # dim(exposure.GWAS) 22 13

# Extract scientific notations using regex 
## Only 2 values to extract for p1: 1e-5, 5e-8
## Only 2 values to extract for p2: 1e-5, 1e-6
library(dplyr, lib.loc = "/software/R/R-3.4.1/lib64/library")

exposure.GWAS <- exposure.GWAS %>%
  dplyr::mutate(clumping.p1.value= stringr::str_extract(fileNameSpecial,"[1|5]e[-][5|8]")
                ,clumping.p2.value= stringr::str_extract(fileNameSpecial,"[1]e[-][5|6]")) # dim(exposure.GWAS) 22 13

# Create outcome file information
outcome.GWAS <- data.frame(filePath=c( filePath.QCed.GWAS.GSCAN
                                      ,filePath.QCed.GWAS.ICC
                                      ,filePath.QCed.GWAS.UKB)
                           ,consortium=c(rep("GSCAN",times=length(filePath.QCed.GWAS.GSCAN))
                                         ,rep("ICC",times=length(filePath.QCed.GWAS.ICC))
                                         ,rep("UKB",times=length(filePath.QCed.GWAS.UKB)))
                           ,trait=c( "AI","CPD","DPW","SC","SI"
                                    ,"CI"
                                    ,"CCPD","caffeine","ESDPW","PYOS","CPD")
                           ,delimiter=c( rep("tab",times=length(filePath.QCed.GWAS.GSCAN))
                                        ,rep("space",times=length(filePath.QCed.GWAS.ICC))
                                        ,rep("tab",times=length(filePath.QCed.GWAS.UKB)))
                           ,colname.for.SNP=c( rep("RSID",times=length(filePath.QCed.GWAS.GSCAN))
                                              ,rep("SNP",times=length(filePath.QCed.GWAS.ICC))
                                              ,rep("SNP",times=length(filePath.QCed.GWAS.UKB)))
                           ,colname.for.beta=c( rep("BETA",times=length(filePath.QCed.GWAS.GSCAN))
                                                ,rep("Effect",times=length(filePath.QCed.GWAS.ICC))
                                                ,rep("BETA",times=length(filePath.QCed.GWAS.UKB)))
                           ,colname.for.SE=c( rep("SE",times=length(filePath.QCed.GWAS.GSCAN))
                                              ,rep("StdErr",times=length(filePath.QCed.GWAS.ICC))
                                              ,rep("SE",times=length(filePath.QCed.GWAS.UKB)))
                           ,colname.for.effect.allele=c( rep("ALT",times=length(filePath.QCed.GWAS.GSCAN))
                                                         ,rep("Allele1",times=length(filePath.QCed.GWAS.ICC))
                                                         ,rep("ALLELE1",times=length(filePath.QCed.GWAS.UKB)))
                           ,colname.for.other.allele=c( rep("REF",times=length(filePath.QCed.GWAS.GSCAN))
                                                        ,rep("Allele2",times=length(filePath.QCed.GWAS.ICC))
                                                        ,rep("ALLELE0",times=length(filePath.QCed.GWAS.UKB)))
                           ,colname.for.P.value=c( rep("PVALUE",times=length(filePath.QCed.GWAS.GSCAN))
                                                   ,rep("P",times=length(filePath.QCed.GWAS.ICC))
                                                   ,rep("PVALUE",times=length(filePath.QCed.GWAS.UKB)))
                           ,stringsAsFactors = FALSE) # dim(outcome.GWAS) 11 10

# Export files
ExportFileTabSeparated( data=exposure.GWAS
                       ,output.file.path = paste0(loc.twoSampleMR.input,"file-info_exposure-clumped-GWASs.tsv"))

ExportFileTabSeparated( data=outcome.GWAS
                        ,output.file.path = paste0(loc.twoSampleMR.input,"file-info_outcome-QCed-GWASs.tsv"))

#file.copy(paste0(locScripts,"MR_step06-03-01_create-iterators-for-running-two-sample-MR-on-HPC.R"),paste0(locScripts,"MR_step08-01-01_create-files-to-loop-through-running-LDSC.R"))
#-------------------------------------------------------------------------------------------------
#------------------This is the end of this file---------------------------------------------------
#-------------------------------------------------------------------------------------------------
