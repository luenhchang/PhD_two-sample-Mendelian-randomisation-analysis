#!/usr/bin/env Rscript

##################################################################################
# Filename: MR_step08-02-01_create-files-to-loop-through-running-LDSC-genetic-correlations.R
# Modified from: MR_step08-01-01_create-files-to-loop-through-running-LDSC-munge.R
# Program author: Chang
# Purpose: Create a tsv file for next step to loop through calculating genetic correlations using LDSC software
# Date created: 20190411
# Note: missing values should be "nan" for LDSC software
# File directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Input paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv")
# Outpu paste0(loc.LDSC.input,"file-info_munged-QCed-GWASs.tsv")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20191001  Exported the 1 files above
# 20190813  Exported the 1 files above
# 20190414  Exported the 1 files above
# 20190411  Exported the 1 files above
#-----------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

#-------------------------------------------
# Folders under lunC working
#-------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
loc.LDSC <- paste0(workingDir,"MR_ICC_GSCAN_201806/LD-score-correlation/")
loc.LDSC.input <- paste0(loc.LDSC,"input/")
loc.LDSC.munged <- paste0(loc.LDSC,"output/munged-GWASs")

#--------------------------------------------------------------------
# Folder locations for output files
#--------------------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

# Import a tsv file containing information about QCed GWAS files
ImportATabSeparatedFile(input.file.path = paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv")
                        , data.name = "QCed.GWASs" ) # dim(QCed.GWASs) 11 17

#------------------------------------------------------------------------------
# Get file paths of Munged QCed GWAS files
#------------------------------------------------------------------------------
munged.files <- data.frame(filePath=list.files(loc.LDSC.munged
                                               ,pattern = glob2rx("^HapMap3-SNPs-found-in_QCed-GWAS-*-*_munged.sumstats.gz$")
                                               ,full.names = TRUE)
                           ,stringsAsFactors = F) # dim(munged.files) 11 1

munged.files$fileName <- with(munged.files, basename(filePath)) # dim(munged.files) 11 2

# Extract consortium and trait names from file paths
munged.files2 <- munged.files %>% tidyr::separate(col=fileName
                                                  ,into=c("fileName_pt1","fileName_pt2","fileName_pt3")
                                                  ,sep="_"
                                                  ,remove=TRUE) %>%
  tidyr::separate(col=fileName_pt2
                  ,into=c("fileName_pt2.1","fileName_pt2.2","consortium","trait")
                  ,sep="-"
                  ,remove=TRUE) %>%
  dplyr::mutate(substance= dplyr::case_when(trait %in% c("AI","CPD","PYOS","SC","SI") ~ "tobacco"
                                     , trait %in% c("DPW","ESDPW") ~ "alcohol"
                                     , trait %in% c("CI") ~ "cannabis"
                                     , trait %in% c("CCPD") ~ "caffeine"
                                     , TRUE ~ as.character(trait))) %>%
  dplyr::select(-one_of(c("fileName_pt1","fileName_pt3","fileName_pt2.1","fileName_pt2.2"))) # dim(munged.files2) 11 4


# Merge munged files with sample size and prevalence of QCed GWAS
munged.files3 <- dplyr::left_join(munged.files2
                                  ,QCed.GWASs
                                  ,by=c("consortium" = "consortium"
                                        ,"trait" = "trait")) %>%
  dplyr::select_(.dots = c("filePath.x","consortium","trait","substance","trait.type","N.continuous","N.cases","N.controls","population.prevalence","trait.definition","data.type")) # dim(munged.files3) 11 11

# Calculate sample prevalence using N.cases and N.controls
munged.files3$sample.prevalence <- with(munged.files3, N.cases/(N.cases+N.controls)) # dim(munged.files3) 11 12

# Export files
ExportFileTabSeparated( data=munged.files3
                        ,missing.values.as = "nan"
                        ,output.file.path = paste0(loc.LDSC.input,"file-info_munged-QCed-GWASs.tsv"))

#file.copy(paste0(locScripts,"MR_step08-02-01_create-files-to-loop-through-running-LDSC-genetic-correlations.R"))

#-------------------------------------------------------------------------------------------------
#------------------This is the end of this file---------------------------------------------------
#-------------------------------------------------------------------------------------------------
