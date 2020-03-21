#!/usr/bin/env Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step06-03-06_run_MR-leave-one-out_Egger-intercept.R
# Modified from : MR_step06-03-05_tabulate-two-sample-MR-analysis-results.R zMR_step06-03_two-sample-MR.R
# Date created  : 20190709
# Purpose       : Run MR leave one out analysis, pleiotropy test (Egger intercept) on selective pairs of exposures and outcomes                 
# Note          : 
#----------------------------------------------------------------------------------------
# Run dependency:     
# Function external:  
# Type  Files
#----------------------------------------------------------------------------------------
# Input paste0(loc.twoSampleMR.harmonised,"harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine.tsv")
# Input paste0(loc.twoSampleMR.harmonised,"harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5_outcome-UKB-caffeine.tsv")

# Outpu paste0(loc.twoSampleMR.tabulated,"leave-one-SNP-out-analysis-results_association-GSCAN-smoking-initiation_on_UKB-caffeine-consumed-per-day.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190823 Exported paste0(loc.twoSampleMR.tabulated,"leave-one-SNP-out-analysis-results_association-GSCAN-smoking-initiation_on_UKB-caffeine-consumed-per-day.tsv")
# 20190821 Converted causal estimates to represent a doubling in odds of binary exposure, exported paste0(loc.twoSampleMR.tabulated,"leave-one-SNP-out-analysis-results_association-GSCAN-smoking-initiation_on_UKB-caffeine-consumed-per-day.tsv")
# 20190819  Exported paste0(loc.twoSampleMR.tabulated,"leave-one-SNP-out-analysis-results_association-GSCAN-smoking-initiation_on_UKB-caffeine-consumed-per-day.tsv")
# 20190709  Ran MR leave one out and MR IVW. Results pasted to emails to Ong JS
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
loc.twoSampleMR.harmonised <- paste0(loc.twoSampleMR,"input/harmonised-data/")
loc.twoSampleMR.tabulated <- paste0(loc.twoSampleMR,"result-tabulated/")

#-----------------------------------------------------------------------------------------
# Load packages, import external functions
#-----------------------------------------------------------------------------------------
library("TwoSampleMR",lib.loc="/software/R/R-3.5.1/lib64/R/library")
#library("dplyr",lib.loc="/software/R/R-3.5.1/lib64/R/library")
library("magrittr",,lib.loc="/software/R/R-3.5.1/lib64/R/library")
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#------------------------------------------------------------------------------
# Import harmonised data in TSV format
#------------------------------------------------------------------------------
columns.keep <- c("exposure","outcome","SNP","b","se","p")

ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.harmonised,"harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine.tsv")
                        ,data.name="harmonised.GSCAN.SI.5e8.UKBECCPD") # dim(harmonised.GSCAN.SI.5e8.UKBECCPD) 7 28

ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.harmonised,"harmonised-data_exposure-clumped-si-noICC-LDWindow-kb-10000-R2-0.01-p1-1e-5-p2-1e-5_outcome-UKB-caffeine.tsv")
                        ,data.name="harmonised.GSCAN.SI.1e5.UKBECCPD") # dim(harmonised.GSCAN.SI.1e5.UKBECCPD) 67 28

# Examine CPD and caffeine 
ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.harmonised,"harmonised-data_exposure-clumped-cpd-noICC-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine.tsv")
                        ,data.name="harmonised.GSCAN.CPD.5e8.UKBECCPD") # dim(harmonised.GSCAN.CPD.5e8.UKBECCPD) 10 28

ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.harmonised,"harmonised-data_exposure-clumped-GWAS-UKB3456-LDWindow-kb-10000-R2-0.01-p1-5e-8-p2-1e-6_outcome-UKB-caffeine.tsv")
                        ,data.name="harmonised.UKB.CPD.5e8.UKBECCPD") # dim(harmonised.UKB.CPD.5e8.UKBECCPD) 3 28

#------------------------------------------------------------------------------
# Perform MR leave-one-out analysis, running MR by removing one SNP at a time, to identify if a single SNP is driving the association.
#------------------------------------------------------------------------------
# Scale up causal estimate of caffeine exposure up using 1 cup change of coffee (75mg) or tea (40 m) consumption
caffeine.per.cup.coffee <- 75
caffeine.per.cup.tea <- 40
caffeine.per.cup <- (caffeine.per.cup.coffee+caffeine.per.cup.tea)/2

# Run leave-one-out on exposure UKB-CCPD and outcome ICC-CI
LOO1 <- TwoSampleMR::mr_leaveoneout(harmonised.GSCAN.SI.5e8.UKBECCPD) %>% 
  dplyr::select_(.dots=columns.keep) %>%
  dplyr::mutate( SNP=as.character(SNP)
                ,var.type.exposure=case_when( exposure %in% c("SC","SI","CI") ~ "binary"
                                              ,exposure %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                ,exposure.substance=case_when( exposure %in% c("AI","CPD","PYOS","SC","SI") ~ "tobacco"
                                               ,exposure %in% c("DPW","ESDPW") ~ "alcohol"
                                               ,exposure %in% c("CI") ~ "cannabis"
                                               ,exposure %in% c("caffeine") ~ "caffeine"
                                               ,TRUE ~ as.character(exposure))
                ,var.type.outcome=case_when( outcome %in% c("SC","SI","CI") ~ "binary"
                                             ,outcome %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                ,outcome.substance=case_when( outcome %in% c("AI","CPD","PYOS","SC","SI") ~ "tobacco"
                                              ,outcome %in% c("DPW","ESDPW") ~ "alcohol"
                                              ,outcome %in% c("CI") ~ "cannabis"
                                              ,outcome %in% c("caffeine") ~ "caffeine"
                                              ,TRUE ~ as.character(outcome))
                ,effect.size.estimate=case_when( 
                   var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b
                  ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ b*caffeine.per.cup
                  ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b)
                  ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(b*caffeine.per.cup)
                  ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*b
                  ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*b)
                                                )
                ,effect.size.LBound=case_when(
                   var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b-1.96*se
                  ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b-1.96*se)
                  ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b-1.96*se)
                  ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b-1.96*se))
                  ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b-1.96*se)
                  ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b-1.96*se))
                                              )
                ,effect.size.UBound=case_when(
                   var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b+1.96*se
                  ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b+1.96*se)
                  ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b+1.96*se)
                  ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b+1.96*se))
                  ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b+1.96*se)
                  ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b+1.96*se))
                                              )
                # Correct for multiple testing. Use Bonferroni correction
                ,sig.thres.BC=0.05/nrow(.)
                ,significant=case_when(p < sig.thres.BC ~ "yes"
                                           , TRUE ~ "no")
                ) # dim(LOO1) 8 15

LOO2 <- TwoSampleMR::mr_leaveoneout(harmonised.GSCAN.SI.1e5.UKBECCPD) %>% 
  dplyr::select_(.dots=columns.keep) %>%
  dplyr::mutate( SNP=as.character(SNP)
                 ,var.type.exposure=case_when( exposure %in% c("SC","SI","CI") ~ "binary"
                                               ,exposure %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                 ,exposure.substance=case_when( exposure %in% c("AI","CPD","PYOS","SC","SI") ~ "tobacco"
                                                ,exposure %in% c("DPW","ESDPW") ~ "alcohol"
                                                ,exposure %in% c("CI") ~ "cannabis"
                                                ,exposure %in% c("caffeine") ~ "caffeine"
                                                ,TRUE ~ as.character(exposure))
                 ,var.type.outcome=case_when( outcome %in% c("SC","SI","CI") ~ "binary"
                                              ,outcome %in% c("AI","CPD","DPW","caffeine","ESDPW","PYOS") ~ "continuous")
                 ,outcome.substance=case_when( outcome %in% c("AI","CPD","PYOS","SC","SI") ~ "tobacco"
                                               ,outcome %in% c("DPW","ESDPW") ~ "alcohol"
                                               ,outcome %in% c("CI") ~ "cannabis"
                                               ,outcome %in% c("caffeine") ~ "caffeine"
                                               ,TRUE ~ as.character(outcome))
                 ,effect.size.estimate=case_when( 
                    var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b
                   ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ b*caffeine.per.cup
                   ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b)
                   ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(b*caffeine.per.cup)
                   ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*b
                   ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*b)
                                                )
                 ,effect.size.LBound=case_when(
                    var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b-1.96*se
                   ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b-1.96*se)
                   ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b-1.96*se)
                   ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b-1.96*se))
                   ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b-1.96*se)
                   ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b-1.96*se))
                                              )
                 ,effect.size.UBound=case_when(
                    var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="continuous" ~ b+1.96*se
                   ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="continuous" ~ caffeine.per.cup*(b+1.96*se)
                   ,var.type.exposure=="continuous" & exposure !="caffeine" & var.type.outcome=="binary" ~ exp(b+1.96*se)
                   ,var.type.exposure=="continuous" & exposure =="caffeine" & var.type.outcome=="binary" ~ exp(caffeine.per.cup*(b+1.96*se))
                   ,var.type.exposure=="binary" & var.type.outcome=="continuous" ~ log(2)*(b+1.96*se)
                   ,var.type.exposure=="binary" & var.type.outcome=="binary" ~ exp(log(2)*(b+1.96*se))
                                              )
                 # Correct for multiple testing. Use Bonferroni correction
                 ,sig.thres.BC=0.05/nrow(.)
                 ,significant=case_when(p < sig.thres.BC ~ "yes"
                                        , TRUE ~ "no")
                 ) # dim(LOO2) 68 15

# Join the two tables
LOO.GSCAN.SI.on.UKB.caffeine <- dplyr::full_join(LOO2,LOO1
                                                 ,by=c("SNP"="SNP"
                                                       ,"exposure"="exposure"
                                                       ,"var.type.exposure"="var.type.exposure"
                                                       ,"outcome"="outcome"
                                                       ,"var.type.outcome"="var.type.outcome")
) # dim(LOO.GSCAN.SI.on.UKB.caffeine) 68 25

# Replace .x and .y colname suffix with clumping p1 value
colnames.new <- gsub(colnames(LOO.GSCAN.SI.on.UKB.caffeine),pattern = "\\.x",replacement = "\\.1e5")
colnames(LOO.GSCAN.SI.on.UKB.caffeine) <- gsub(colnames.new, pattern = "\\.y",replacement = "\\.5e8")

# Export data
ExportFileTabSeparated(data = LOO.GSCAN.SI.on.UKB.caffeine
                        ,missing.values.as = ""
                        ,output.file.path = paste0(loc.twoSampleMR.tabulated,"leave-one-SNP-out-analysis-results_association-GSCAN-smoking-initiation_on_UKB-caffeine-consumed-per-day.tsv"))

#file.copy(from=paste0(locScripts,"MR_step06-03-06_run_MR-leave-one-out_MR-IVW.R"),to=paste0(locScripts,"MR_step06-03-07_run_heterogeneity-test.R"))

#----------------------------------------------------------------------------#
# Run horizontal pleiotropy test
#----------------------------------------------------------------------------#
TwoSampleMR::mr_pleiotropy_test(dat=harmonised.GSCAN.SI.5e8.UKBECCPD)
#   id.exposure id.outcome  outcome exposure egger_intercept      se      pval
# 1      70JK0x     lY7d1Y caffeine       SI      -0.5113987 2.49741 0.8458255

TwoSampleMR::mr_pleiotropy_test(dat=harmonised.GSCAN.SI.1e5.UKBECCPD)
#   id.exposure id.outcome  outcome exposure egger_intercept        se      pval
# 1      fB9Dlg     OC9Xxe caffeine       SI       0.1525081 0.3784044 0.6882492

#----------------------------------------------------------------------------#
# CPD and caffeine consumption
#----------------------------------------------------------------------------#
unique(harmonised.GSCAN.CPD.5e8.UKBECCPD$SNP)
#[1] "rs193122850" "rs2118359"   "rs215600"    "rs2273500"   "rs56113850"  "rs62012628"  "rs7017612"   "rs73467256"  "rs7872903" "rs9788721" 
unique(harmonised.UKB.CPD.5e8.UKBECCPD$SNP)
# "rs112878080" "rs56113850"  "rs774020296"

# Subset exposure SNP rs16969968 from GSCAN GWAS files in Unix 
# cd /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/noICC_results/QC3_remove_ambiguousSNPs_indel
# cat <(head -1 cpd_noICC.ambiguousSNPRemoved) <(grep rs16969968 cpd_noICC.ambiguousSNPRemoved) > cpd_noICC.ambiguousSNPRemoved_rs16969968
file.path.rs16969968.GSCAN.CPD <- "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/noICC_results/QC3_remove_ambiguousSNPs_indel/cpd_noICC.ambiguousSNPRemoved_rs16969968"

# Subset exposure SNP rs16969968 from UKB GWAS files in Unix 
# cd /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel
# cat <(head -1 QCed-GWAS-UKB3456_headed) <(grep rs16969968 QCed-GWAS-UKB3456_headed) > QCed-GWAS-UKB3456_headed_rs16969968
file.path.rs16969968.UKB.CPD <- "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB3456_headed_rs16969968"

# Subset outcome SNP rs16969968 from UKB GWAS on caffeine consumption in Unix 
# cd /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC3_remove_ambiguousSNPs_indel
# cat <(head -1 QCed-GWAS-UKB-caffeine-consumed-per-day_headed) <(grep rs16969968 QCed-GWAS-UKB-caffeine-consumed-per-day_headed) > QCed-GWAS-UKB-caffeine-consumed-per-day_headed_rs16969968
file.path.rs16969968.UKB.caffeine <- "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-caffeine-consumed-per-day_headed_rs16969968"

ImportATabSeparatedFile(input.file.path = file.path.rs16969968.GSCAN.CPD,data.name = "rs16969968.GSCAN.CPD")
ImportATabSeparatedFile(input.file.path = file.path.rs16969968.UKB.CPD,data.name = "rs16969968.UKB.CPD")
ImportATabSeparatedFile(input.file.path = file.path.rs16969968.UKB.caffeine,data.name = "rs16969968.UKB.caffeine")

# Format data as exposure data 
rs16969968.GSCAN.CPD.exposure <- TwoSampleMR::format_data( dat=rs16969968.GSCAN.CPD
                                                          ,type="exposure"
                                                          ,header = TRUE
                                                          ,snp_col = "RSID"
                                                          ,beta_col = "BETA"
                                                          ,se_col = "SE"
                                                          ,effect_allele_col = "ALT"
                                                          ,other_allele_col = "REF"
                                                          ,pval_col = "PVALUE")

rs16969968.UKB.CPD.exposure <- TwoSampleMR::format_data( dat=rs16969968.UKB.CPD
                                                           ,type="exposure"
                                                           ,header = TRUE
                                                           ,snp_col = "SNP"
                                                           ,beta_col = "BETA"
                                                           ,se_col = "SE"
                                                           ,effect_allele_col = "ALLELE1"
                                                           ,other_allele_col = "ALLELE0"
                                                           ,pval_col = "PVALUE")

rs16969968.UKB.caffeine.outcome <- TwoSampleMR::format_data( dat=rs16969968.UKB.caffeine
                                                             ,type="outcome"
                                                             ,header = TRUE
                                                             ,snp_col = "SNP"
                                                             ,beta_col = "BETA"
                                                             ,se_col = "SE"
                                                             ,effect_allele_col = "ALLELE1"
                                                             ,other_allele_col = "ALLELE0"
                                                             ,pval_col = "PVALUE")
# Inner join exposure and outcome
harmo.data.GSCAN.CPD.UKB.caffeine <- TwoSampleMR::harmonise_data( exposure_dat=rs16969968.GSCAN.CPD.exposure
                                                                  ,outcome_dat = rs16969968.UKB.caffeine.outcome)

harmo.data.UKB.CPD.UKB.caffeine <- TwoSampleMR::harmonise_data( exposure_dat=rs16969968.UKB.CPD.exposure
                                                                  ,outcome_dat = rs16969968.UKB.caffeine.outcome)

# Run single SNP MR
TwoSampleMR::mr(harmo.data.GSCAN.CPD.UKB.caffeine)
#   id.exposure id.outcome outcome exposure     method nsnp        b      se     pval
# 1      6Tq93o     tEsqqK outcome exposure Wald ratio    1 4.582536 6.84695 0.503316

TwoSampleMR::mr(harmo.data.UKB.CPD.UKB.caffeine)
#   id.exposure id.outcome outcome exposure     method nsnp        b        se     pval
# 1      F2EEag     tEsqqK outcome exposure Wald ratio    1 0.439979 0.6573901 0.503316

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

