#!/usr/bin/Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step09-01_make-table1-data-sources-sample-sizes-types-of-analysis-performed-on-substance-use-traits.R
# Modified from : PRS_UKB_201711_step22-03_heatmap-genetic-correlations.R
# Date created  : 20190725
# Purpose       : Make table 1 using Table 1 from the manuscript
# Note: 
#----------------------------------------------------------------------------------------
# Run dependency: 
# MR_step08-03_parse-tabulate_LDSC-SNP-heritability_LDSC-genetic-correlations.R
# MR_step08-02-01_create-files-to-loop-through-running-LDSC-genetic-correlations.R
# MR_step08-01-01_create-files-to-loop-through-running-LDSC-munge.R

# Function external: ExportFileTabSeparated()

# Type  Files
#----------------------------------------------------------------------------------------------
# Input D:\googleDrive\GoogleDriveAsMaster\Chang_manuscript04_MR-cannabis-licitSubstanceUse\Addiction\manuscript4-main-text[2]citation-style-Vancouver-modified.docx
# Input paste0(loc.LDSC.tabulated,"LDSC-SNP-heritability_sample-size-prevalence.tsv")
# Input paste0(loc.obs.assoc,"sample-sizes_phenotypic-associations.tsv")

# Outpu paste0(output.folder,"manu4_data-sources_cohorts_sample-sizes_analysis-performed.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20191010, 20191001, 20190904, 20190815, 20190801, 20190725  
# Exported manu4_data-sources_cohorts_sample-sizes_analysis-performed.tsv
#----------------------------------------------------------------------------------------

#---------------------------------------------
# Folder locations under my home directory
#---------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
paste0(loc.obs.assoc,"sample-sizes_phenotypic-associations.tsv")

#---------------------------------------------
# Folder locations under my working directory
#---------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
loc.LDSC <- paste0(workingDir,"MR_ICC_GSCAN_201806/LD-score-correlation/")
loc.LDSC.input <- paste0(loc.LDSC,"input/")
loc.LDSC.tabulated <- paste0(loc.LDSC,"output/result-tabulated/");
loc.obs.assoc <- paste0(workingDir,"MR_ICC_GSCAN_201806/observational-associations/")
output.folder <- paste0(workingDir,"MR_ICC_GSCAN_201806/study-information/")

#dir.create(output.folder)

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#---------------------------------------------
# Import data
#---------------------------------------------
# Sample sizes from phenotypic associations
# obs.sample.sizs <- ImportATabSeparatedFile(input.file.path = paste0(loc.obs.assoc,"sample-sizes_phenotypic-associations.tsv")
#                                            ,data.name = "obs.sample.sizs") %>%
#   dplyr::mutate(consortium="UKB"
#                 ,dep.var.label=dplyr::recode(dep.var.label,`ECCPD`="caffeine")
#                 ,data.type="Phenotype")
  # dim(obs.sample.sizs) 22 8

# SNP heritability estimates and sample sizes of GWASs
SNP.h2.sample.sizes <- ImportATabSeparatedFile(input.file.path = paste0(loc.LDSC.tabulated,"LDSC-SNP-heritability_sample-size-prevalence.tsv")
                                               ,data.name = "SNP.h2.sample.sizes") %>%
  #dplyr::mutate(substance=dplyr::recode(substance,`tobacco`="nicotine")) %>%
  #dplyr::mutate(analysis.performed=dplyr::case_when( consortium=="UKB" & trait %in% c("PYOS","ESDPW","caffeine") ~ "OA, GA"
  #                                            ,TRUE ~ "GA")) %>%
  dplyr::rename(sample.size.GWAS=sample.size) %>%
  dplyr::select_(.dots = c("consortium","trait","trait.definition","data.type","trait.type","substance","sample.size.GWAS")) # dim(SNP.h2.sample.sizes) 10 7 ,"analysis.performed

#---------------------------------------------------------------------------------------------------
# Create a data.frame using Table 1 from manuscript4 main text
#---------------------------------------------------------------------------------------------------
UKB.cannabis.initiation.info <- data.frame( consortium="UKB"
                                           ,trait="CI"
                                           ,trait.definition="Ever versus never taken cannabis (UKB Data-field: 20453)"
                                           ,data.type="Phenotype"
                                           ,trait.type="binary"
                                           ,substance="cannabis"
                                           #,analysis.performed="OA"
                                           ,stringsAsFactors = F) # dim(UKB.cannabis.initiation.info) 1 7 #,sample.size=153726

UKB.smoking.initiation.info <- data.frame( consortium="UKB"
                                            ,trait="SI"
                                           ,trait.definition="Ever versus never smoked (UKB Data-field: 20160)"
                                           ,data.type="Phenotype"
                                            ,trait.type="binary"
                                            ,substance="nicotine"
                                            #,analysis.performed="OA"
                                            ,stringsAsFactors = F) # dim(UKB.smoking.initiation.info) 1 7
# Combine the two data sets
all <- dplyr::bind_rows(SNP.h2.sample.sizes
                        ,UKB.cannabis.initiation.info
                        ,UKB.smoking.initiation.info ) # dim(all) 12 7

# Combine phenotypic associations and GWASs
# obs.GWAS <- dplyr::full_join( x=all
#                               ,y=obs.sample.sizs[,c("consortium","dep.var.label","data.type","numb.obs.analysed")]
#                               ,by=c("consortium"="consortium","trait"="dep.var.label")) %>%
#   dplyr::rename(data.type.GWAS=data.type.x, data.type.pheno=data.type.y) # dim(obs.GWAS) 12 10

# Export data as a TSV
ExportFileTabSeparated(data = all
                       , missing.values.as = ""
                       , output.file.path = paste0(output.folder,"manu4_data-sources_cohorts_sample-sizes_analysis-performed.tsv"))
