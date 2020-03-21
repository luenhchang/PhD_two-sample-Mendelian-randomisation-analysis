#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step06-03-05_tabulate-two-sample-MR-analysis-results.R
# Modified from : MR_step06-03-02_jobScript-harmonise-exposure-outcome.R
# Date created  : 20190409
# Purpose       : Combine MR analysis results on multiple pairs of exposure and outcome as one file. The output file has results from MR-Egger, weighted median, weighted mode and inverse variance weighted
#                 Convert beta from MR estimators to interpretable effect size and confidence interval. Caffeine consumed per day is measured in mg. This is scaled up to (75+40)/2 to represent one cup change of exposure
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency:     
# Function external:  
# Type  Files
#----------------------------------------------------------------------------------------
# Input Sys.glob(paste0(loc.twoSampleMR.output,"MR-analysis_exposure-*_outcome-*.tsv")) # 180 files

# Outpu paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# Outpu paste0(loc.twoSampleMR.tabulated,"MR-analysis-results.tsv")
# Outpu paste0(loc.phenotypic.association,"exposure-outcome-pairs-used-in-MR-analyses.tsv")

#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20191024 Exported paste0(loc.phenotypic.association,"exposure-outcome-pairs-used-in-MR-analyses.tsv")
# 20191023  Removed exposure ESDPW, caffeine from 4th status "MR to be stratified by smoking, drinking or caffeine use status in outcome GWAS." Exported  MR-analysis-results_all-trait-pairs.tsv 
# 20191010  Added exposure smoking cessation to 4th status "MR to be stratified by smoking, drinking or caffeine use status in outcome GWAS." Exported  MR-analysis-results_all-trait-pairs.tsv 
# 20191002  Mark exposures that were measured in smokers (e.g. CPD, PYOS), drinkers (DPW, ESDPW) or caffeine users (CCPD, ECCPD). MR can be done for these only if there is individual level information to stratify outcome GWAS by smoking, drinking or caffeine use status. Exporte MR-analysis-results_all-trait-pairs.tsv
# 20190926  Exporte MR-analysis-results_all-trait-pairs.tsv
# 20190926  Categorised exposure-outcome pairs as (1) from overlapped samples, (2) are logically impossible sequences (e.g. exposure=smoking cessation & outcome= cigarettes per day and more), or (3) identical/similar traits, (4) reasonable two-sample MR. Use (4) as results  
# 20190904 Changed tobacco to nicotine in substance. Exported paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# 20190822 Converted causal estimates from caffeine per day exposure by multiplying beta by caffeine per cup (75+40)/2. 
# 20190822 Exported paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# 20190821 Converted causal estimates to represent a doubling in odds of binary exposure, exported paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# 20190816 Exported paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# 20190813 Exported paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# 20190806 Exported paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")
# 20190420 Exported the 1 file above
# 20190408 Exported the 1 file above
#----------------------------------------------------------------------------------------

#---------------------------------------------------------------------
# Folder locations under my home directory
#---------------------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

#---------------------------------------------------------------------
# Folder locations under my working directory
#---------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
loc.twoSampleMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/two-sample-MR")
loc.twoSampleMR.output <- paste0(loc.twoSampleMR,"/output/")
loc.twoSampleMR.tabulated <- paste0(loc.twoSampleMR,"/result-tabulated/")
loc.twoSampleMR.report <- paste0(loc.twoSampleMR,"/output_MR-reports")
loc.twoSampleMR.report.figure <- paste0(loc.twoSampleMR,"/output_MR-reports/figure")
loc.twoSampleMR.pleiotropy <- paste0(loc.twoSampleMR,"/output_horizontal-pleiotropy")
loc.phenotypic.association <- paste0(workingDir,"MR_ICC_GSCAN_201806/observational-associations/")

#dir.create(loc.twoSampleMR.tabulated)

#----------------------------------------------------------
# Import functions
#----------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#---------------------------------------------------------------------
# Create a data.frame containing info about individual MR result files
#---------------------------------------------------------------------
# Get a vector of file paths of all MR result files in the input folder
filePath.all.MR.files <- Sys.glob(paste0(loc.twoSampleMR.output,"MR-analysis_exposure-*.tsv")) # length(filePath.all.MR.files) 220
fileName.all.MR.files <- basename(filePath.all.MR.files) # length(fileName.all.MR.files) 220

all.MR <- data.frame(filePath=filePath.all.MR.files
                     ,fileName=fileName.all.MR.files
                     ,stringsAsFactors = F) # dim(all.MR) 220 2

# Split the fileName column into columns of exposure consortium, trait and outcome consortium and trait
all.MR <- all.MR %>% tidyr::separate(col=fileName
                                     ,into=c("fileName.prefix","fileName.exposure","fileName.outcome")
                                     ,sep= "_"
                                     ,remove=TRUE) %>% 
  tidyr::separate(col=fileName.exposure
                  ,into=c("exposure.prefix","exposure.consortium","exposure.trait","exposure.p1.part1","exposure.p1.part2","exposure.p2.part1","exposure.p2.part2")
                  ,sep="-"
                  ,remove=TRUE) %>%
  tidyr::separate(col=fileName.outcome
                  ,into=c("outcome.prefix","outcome.consortium","outcome.trait")) %>%
  dplyr::select(-one_of(c("fileName.prefix","exposure.prefix","outcome.prefix"))) # dim(all.MR) 220 9

# Categorise MR exposure-outcome trait pairs as (1) overlapped samples, (2) identical traits, (3) impossible temporal sequence, (4) continuous exposures measured in smokers, drinkers, or caffeine users, which MR is not possibly done by stratifying outcome GWAS by smoking, drinking or caffeine use status.
all.MR2 <- all.MR %>% 
  dplyr::mutate(status=dplyr::case_when(  exposure.consortium==outcome.consortium ~ "overlapped samples"
                                         ,(exposure.trait==outcome.trait)|(exposure.trait=="DPW" & outcome.trait=="ESDPW")|(exposure.trait=="ESDPW" & outcome.trait=="DPW") ~ "identical traits"
                                         ,(exposure.trait=="SC" & outcome.trait %in% c("AI","SI","CPD"))|(exposure.trait=="CPD" & outcome.trait %in% c("AI","SI","PYOS"))|(exposure.trait=="PYOS" & outcome.trait %in% c("CPD","SI","AI")) ~ "logically impossible sequence"
                                         ,(exposure.trait %in% c("AI","SC","CPD","DPW","CCPD","PYOS")) ~ "MR to be stratified by smoking, drinking or caffeine use status in outcome GWAS"
                                         ,outcome.trait=="CCPD" ~ "Traits not in use"
                                         ,TRUE ~ "reasonable two-sample MR"
                                         ))
table(all.MR2$status)

all.MR2.small <- all.MR2 %>%
  dplyr::filter(status=="reasonable two-sample MR" & outcome.trait!="DPW") %>%
  dplyr::select_(.dots=c("outcome.consortium","outcome.trait","exposure.consortium","exposure.trait")) %>%
  dplyr::arrange(outcome.consortium,outcome.trait,exposure.consortium,exposure.trait) %>%
  dplyr::distinct()
  # dim(all.MR2.small) 23 4

ExportFileTabSeparated(data=all.MR2.small
                       , output.file.path = paste0(loc.phenotypic.association,"exposure-outcome-pairs-used-in-MR-analyses.tsv"))


#-------------------------------------------------------------------------------------------------
# Combine all MR result files as a single file using information from the data.frame created above
#-------------------------------------------------------------------------------------------------
# Read and combine multiple MR result files as a file
append_MR_files <- data.frame()

for (line in 1:nrow(all.MR2)){
  filePath <- all.MR2[line,"filePath"] 
  exposure.consortium <- all.MR2[line,"exposure.consortium"]
  exposure.trait <- all.MR2[line,"exposure.trait"]
  outcome.consortium <- all.MR2[line,"outcome.consortium"]
  outcome.trait <- all.MR2[line,"outcome.trait"]
  status <- all.MR2[line,"status"]
  ImportATabSeparatedFile(input.file.path = filePath,data.name = "temp") # dim(temp) 5 13
  temp$exposure.outcome.status <- status
  append_MR_files <- rbind(append_MR_files,temp)
}

# dim(append_MR_files) 1020  14

# Scale up causal estimate of caffeine exposure up using 1 cup change of coffee (75mg) or tea (40 m) consumption
caffeine.per.cup.coffee <- 75
caffeine.per.cup.tea <- 40
caffeine.per.cup <- (caffeine.per.cup.coffee+caffeine.per.cup.tea)/2

# Add abbreviated MR methods, type of exposure, outcome substance
append_MR_files_selective_methods <- append_MR_files %>% 
  # Subset data from 4 MR methods, and exclude CCPD from both exposure and outcome traits
  dplyr::filter(method %in% c("MR Egger","Weighted median","Inverse variance weighted","Weighted mode")) %>%
  dplyr::filter(!(exposure=="CCPD" | outcome=="CCPD")) %>%
  # Add grouping variables
  dplyr::mutate(var.type.exposure=dplyr::case_when( exposure %in% c("SC","SI","CI") ~ "binary"
                                        ,exposure %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                ,exposure.substance=dplyr::case_when( exposure %in% c("AI","CPD","PYOS","SC","SI") ~ "nicotine"
                                               ,exposure %in% c("DPW","ESDPW") ~ "alcohol"
                                               ,exposure %in% c("CI") ~ "cannabis"
                                               ,exposure %in% c("caffeine") ~ "caffeine"
                                               ,TRUE ~ as.character(exposure))
                ,var.type.outcome=dplyr::case_when( outcome %in% c("SC","SI","CI") ~ "binary"
                                        ,outcome %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                ,outcome.substance=dplyr::case_when( outcome %in% c("AI","CPD","PYOS","SC","SI") ~ "nicotine"
                                              ,outcome %in% c("DPW","ESDPW") ~ "alcohol"
                                              ,outcome %in% c("CI") ~ "cannabis"
                                              ,outcome %in% c("caffeine") ~ "caffeine"
                                              ,TRUE ~ as.character(outcome))
                ,method.abbreviated= dplyr::case_when(method=="MR Egger" ~ "Egger"
                                              ,method=="Weighted median" ~ "W Median"
                                              ,method=="Inverse variance weighted" ~ "IVW"
                                              ,method=="Weighted mode" ~ "W Mode"
                                              ,TRUE ~ as.character(method))) %>%
  dplyr::mutate(effect.size.estimate=dplyr::case_when( 
     var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b
    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ b*caffeine.per.cup
    ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b)
    ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(b*caffeine.per.cup)
    ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*b
    ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*b)
                                                )
      ,effect.size.LBound=dplyr::case_when(
         var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b-1.96*se
        ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b-1.96*se)
        ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b-1.96*se)
        ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b-1.96*se))
        ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b-1.96*se)
        ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b-1.96*se))
                                  )
      ,effect.size.UBound=dplyr::case_when(
         var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b+1.96*se
        ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b+1.96*se)
        ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b+1.96*se)
        ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b+1.96*se))
        ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b+1.96*se)
        ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b+1.96*se))
                                    )) # dim(append_MR_files_selective_methods) 657 22

# Reorder columns
#columns.in.this.order <- c("id.exposure","id.outcome","exposure.consortium","exposure.substance","exposure","var.type.exposure","exposure.clumping.p1","exposure.clumping.p2","outcome.consortium","outcome.substance","outcome","var.type.outcome","nsnp","method","method.abbreviated","b","se","pval","effect.size.estimate","effect.size.LBound","effect.size.UBound")

append_MR_files_selective_methods <- append_MR_files_selective_methods %>% 
  dplyr::select(id.exposure,id.outcome,exposure.consortium,exposure.substance,exposure,var.type.exposure,exposure.clumping.p1,exposure.clumping.p2,outcome.consortium,outcome.substance,outcome,var.type.outcome,exposure.outcome.status,nsnp,method,method.abbreviated,b,se,pval,effect.size.estimate,effect.size.LBound,effect.size.UBound) 
# dim(append_MR_files_selective_methods) 657 22

# Export data
ExportFileTabSeparated(data = append_MR_files_selective_methods
                       ,missing.values.as = "NA"
                       ,output.file.path = paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv"))

setwd(locScripts)
#file.copy("MR_step06-03-05_tabulate-two-sample-MR-analysis-results.R","MR_step06-03-06_run-MR-leave-one-out-analysis.R")

#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#