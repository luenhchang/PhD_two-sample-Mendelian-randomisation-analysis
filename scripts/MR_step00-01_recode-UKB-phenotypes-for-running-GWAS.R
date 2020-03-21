# ---------------------------------------------------------------------------------------
# File path     : /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step00-01_recode-UKB-phenotypes-for-running-GWAS.R
# Modified from : /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step00-00_recode_phenotype_for_GWAS.R 
# Date created  : 20190809
# Programmer:     Chang
# Run dependency: D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_phenoUtility\phenoUtility.bat
# To do before this step: (1) Extract phenotype data using Jiyuan's phenoUtility, a GUI software. See D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_phenoUtility\README.txt, (2) Copy output phenotype files to ${dir_ukbPheno}

# Purpose       : (1) find comparable ukb phenotypes to GSCAN's five phenotypes- cigarettes per day (CPD), smoking initiation (SI), smoking cessation (SC), Age at which an individual started smoking regularly (AI), Drinks per week in individuals who are active drinkers (DPW), (2) If a UKB phenotype is continuous, use the average of 3 instances as the phenotype to run GWAS; if it is binary, use instances with consistent responses:

# Note: (1) The columns of the exported phenotype files must in this order: "FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","phenotypeColumn" (2) In plink, Case/control phenotypes are expected to be encoded as 1=unaffected (control), 2=affected (case); 0 is accepted as an alternate missing value encoding. 

#------------------------------------------------------------------------------------------------------
# Input paste0(dir_ukbPheno,"/ukb2907_everStoppedSmokingFor6+Months/ukb2907.phenoUtility")
# Input paste0(dir_ukbPheno,"/ukb3456_numCigareDaily")
# Input paste0(dir_ukbPheno,"/ukb20160_everSmoked/ukb20160.phenoUtility")
# Input paste0(dir_ukb20453,"/ukb20453.phenoUtility")
# Input paste0(lab_stuartma_pheno_alcohol,"total_alcohol_unitsweekly.combined.pheno")
# Input paste0(lab.stuartma.pheno.caffeine,"UKB-data-field-number.phenoUtility")

# Outpu paste0(dir_ukbPheno,"/ukb2907_everStoppedSmokingFor6+Months/ukb2907.phenoUtility.ever_stoppedSmoking")
# Outpu paste0(dir_ukbPheno,"/ukb20160_everSmoked/ukb20160.phenoUtility.recoded")
# Outpu paste0(dir_ukb20453,"/ukb20453.phenoUtility.recoded")
# Outpu paste0(dir_ukb3456,"/ukb3456_IID_NA_in_20453")
# Outpu paste0(outputFolderPath_ukb3456_NA20453,"/ukb3456_IID_NA_in_20453")
# Outpu paste0(outputFolderPath_ESDPW_NA20453,"/phenotype")
# Outpu paste0(outputFolderPath_ukb20161_NA20453,"/phenotype")
# Outpu paste0(outputFolderPath_ukb20161_NA20453,"/phenotype")
# Outpu paste0(lab.stuartma.pheno.caffeine,"daily.caffeine.consumption.thru.coffee.tea")
# Outpu paste0(output.folder.path.UKB.caffeine.NA20453,"phenotype")
# Outpu paste0(dir_ukb20453,file_ukb20453,".recoded")

#------------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20191024    Exported /mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded
# 20190810    Exported paste0(output.folder.path.UKB.caffeine.NA20453,"phenotype")
# 20180820    Exported paste0(outputFolderPath_ESDPW_NA20453,"/phenotype")
# 20180730    Exported paste0(outputFolderPath_ukb3456_NA20453,"/ukb3456_IID_NA_in_20453")
# 20180201    Exported paste0(dir_ukbPheno,"/ukb2907_everStoppedSmokingFor6+Months/ukb2907.phenoUtility.ever_stoppedSmoking")
# 20180201    Checked frequency of paste0(dir_ukbPheno,"/ukb20160_everSmoked/ukb20160.phenoUtility"). Column ever_smoked is taken from X20160.0.0
#----------------------------------------------------------------------------------------
# GSCAN   UKB   type    Phenotype 
# -----------------------------------------------------------------------
# CPD     3456  conti   Number of cigarettes currently smoked daily (current smokers)
# SI      20160	binary  Ever smoked	
# SC      20116 binary  Smoking status
# AI      3436  conti   Age started smoking in current smokers
# DPW     1559  conti   Number of standard drinks per week   
# ------------------------------------------------------------------------

#---------------------------------------------
# Folder paths under my home
#---------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/"
dir_ukbPheno <- paste0(homeDir,"data/UKBiobank_phenotype")
loc.R.functions <- paste0(homeDir,"scripts/RFunctions/")

#---------------------------------------------
# Folder paths under Stuart MC's lab
#---------------------------------------------
lab_stuartmaDir <- "/reference/data/UKBB_500k/versions/lab_stuartma/"
lab_stuartma_pheno <- paste0(lab_stuartmaDir,"pheno/")
lab_stuartma_pheno_smoking <- paste0(lab_stuartma_pheno,"smoking/")
lab_stuartma_pheno_alcohol <- paste0(lab_stuartma_pheno,"alcohol/")
lab.stuartma.pheno.coffee <- paste0(lab_stuartma_pheno,"coffee/")
lab.stuartma.pheno.caffeine <- paste0(lab.stuartma.pheno.coffee,"ukb_1488-teaIntake_1498-coffeeIntake_1508-coffeeType/")

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
lab.stuartma.gwas <- paste0(lab_stuartmaDir,"draft_gwas/");

#---------------------------------------------
# Import external functions
#---------------------------------------------
source(paste0(loc.R.functions,"RFunction_import_export_single_file.R"))

########################################################################################
#--------------------------ukb 20453 Ever taken cannabis ------------------------------#
# http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20453
########################################################################################

dir_ukb20453 <- paste0(dir_ukbPheno,"/ukb20453_everTakenCannabis/")
file_ukb20453 <- "ukb20453.phenoUtility"

# Import phenotype 20453
ImportASpaceSeparatedFile(input.file.path = paste0(dir_ukb20453,file_ukb20453) # dim(pheno_ukb20453) 487409      9
                          ,data.name = "pheno_ukb20453") # dim(pheno_ukb20453) 487409      9

# Check  data values. Only 1 instance is found, not 3
as.data.frame(table(pheno_ukb20453$X20453.0.0, exclude = FALSE))

# Data coding in this phenotype: Data-Coding 526
#------------------------------------------------------------------------------------
#   Var1   Freq Def                       X20453_0_0_recoded        X20453_0_0_recoded_plink
#------------------------------------------------------------------------------------
# 1 -818    225 Prefer not to answer      -818 (uncertain)
# 2    0 119694 No                        0 (never used cannabis)   1 code for controls by plink
# 3    1  14470 Yes, 1-2 times            1 (ever used cannabis)    2 code for cases by plink
# 4    2   8412 Yes, 3-10 times           1 (ever used cannabis)    2  
# 5    3   6824 Yes, 11-100 times         1 (ever used cannabis)    2
# 6    4   4101 Yes, more than 100 times  1 (ever used cannabis)    2
# 7 <NA> 333683 Missing values            NA (no cannabis initiation data)  NA missing by default in plink
#------------------------------------------------------------------------------------

# Categorise X20453.0.0 as: never used cannabis (0), ever used cannabis (1), or no data available (NA)
pheno_ukb20453 <- pheno_ukb20453 %>% 
  # Exclude values: -818 Prefer not to answer
  dplyr::filter(X20453.0.0 %in% c(0,1,2,3,4,NA)) %>%
  dplyr::mutate(X20453_0_0_recoded=dplyr::case_when(X20453.0.0==0 ~ 0
                                             ,X20453.0.0 %in% c(1,2,3,4) ~ 1
                                             ,TRUE ~ as.numeric(NA))) %>% 
    #Recode 0 as 1 for controls, 1 as 2 for cases, and NA as missing as expected by plink
  dplyr::mutate(X20453_0_0_recoded_plink=dplyr::case_when( X20453_0_0_recoded==0 ~ 1
                                                          ,X20453_0_0_recoded==1 ~ 2
                                                          ,TRUE ~ as.numeric(NA))) # dim(pheno_ukb20453) 487184     11
table(pheno_ukb20453$X20453_0_0_recoded, exclude = FALSE)
#      0      1   <NA> 
# 119694  33807 333683

table(pheno_ukb20453$X20453_0_0_recoded_plink, exclude = FALSE)
#      1      2   <NA> 
# 119694  33807 333683 

# Export recoded phenotype back to the pheno folder
ExportFileSpaceSeparated(data=pheno_ukb20453
                         ,missing.values.as="NA"
                         ,output.file.path = paste0(dir_ukb20453,file_ukb20453,".recoded"))

# Subset participants with no cannabis initiation data available. This group will be used for subsetting substance use data from those who have no cannabis data 
pheno_ukb20453_NA <- pheno_ukb20453 %>% 
  dplyr::filter(is.na(X20453_0_0_recoded)) %>%
  dplyr::select_(.dots=c("FID","IID","X20453_0_0_recoded")) # dim(pheno_ukb20453_NA) 333683 3


########################################################################################
#-----Calculate daily caffeine consumption from three phenotypes:----------------------#
#------- 1488 Tea intake (How many cups of tea do you drink each DAY?)
#------- 1498 Coffee intake
#------- 1508 Coffee type
########################################################################################
# Import a phenotype data file
ImportASpaceSeparatedFile(input.file.path = paste0(lab.stuartma.pheno.caffeine,"UKB-data-field-number.phenoUtility")
                          ,data.name = "ukb.pheno.1488.1498.1508") # dim(ukb.pheno.1488.1498.1508) 487409     33

# Check distributions of 3 instances for tea intake
hist(ukb.pheno.1488.1498.1508$X1488.0.0)
hist(ukb.pheno.1488.1498.1508$X1488.1.0)
hist(ukb.pheno.1488.1498.1508$X1488.2.0)
# Coding 100373 defines 3 special values:
# -10 represents "Less than one"
# -1 represents "Do not know"
# -3 represents "Prefer not to answer"

# Check distributions of 3 instances for coffee intake
hist(ukb.pheno.1488.1498.1508$X1498.0.0)
hist(ukb.pheno.1488.1498.1508$X1498.1.0)
hist(ukb.pheno.1488.1498.1508$X1498.2.0)
# Coding 100373 defines 3 special values:
# -10 represents "Less than one"
# -1 represents "Do not know"
# -3 represents "Prefer not to answer"

#-----------------------------------------------------
# Check distributions of 3 instances for 1508 coffee type
#-----------------------------------------------------

# Coding	Meaning                                       Recoding in X1508.0.0.recoded,X1508.1.0.recoded,X1508.2.0.recoded   
#---------------------------------------------------------------------------------------------------------
# 1	      Decaffeinated coffee (any type)               0 (no caffeine consumption thru decaffeinated coffee)
# 2	      Instant coffee                                1 (consume caffeine thru regular coffee)
# 3	      Ground coffee (include espresso, filter etc)  1 (consume caffeine thru regular coffee)
# 4	      Other type of coffee                          1 (consume caffeine thru regular coffee)
# -1	    Do not know                                   NA (no data available)
# -3	    Prefer not to answer                          NA (no data available)
#--------------------------------------------------------------------------------------------------------

table(ukb.pheno.1488.1498.1508$X1508.0.0,exclude = FALSE)
#  -3     -1      1      2      3      4   <NA> 
# 432   1326  72665 210204  86076   6952 109754
table(ukb.pheno.1488.1498.1508$X1508.1.0, exclude = FALSE)
# -3     -1      1      2      3      4   <NA> 
# 11     24   3318   8413   4418    253 470972 
table(ukb.pheno.1488.1498.1508$X1508.2.0,exclude = FALSE)
# -3     -1      1      2      3      4   <NA> 
# 5     12   1910   5122   3229    185 476946 

# Create a function that works like rowMeans() while dealing with NA in dplyr
# URL: https://stackoverflow.com/questions/33401788/dplyr-using-mutate-like-rowmeans/35553114
my_rowmeans = function(..., na.rm=TRUE){
  x = 
    if (na.rm) lapply(list(...), function(x) replace(x, is.na(x), as(0, class(x)))) 
  else       list(...)
  
  d = Reduce(function(x,y) x+!is.na(y), list(...), init=0)
  
  Reduce(`+`, x)/d
} 

# Recode data values from 3 instances, take a mean and calculate daily caffeine consumption
ukb.pheno.1488.1498.1508 <- ukb.pheno.1488.1498.1508 %>% 
  # Recode values <=0 as NA, else as original for the 3 instances
  dplyr::mutate( X1488.0.0.recoded=dplyr::case_when(X1488.0.0 > 0 ~ as.numeric(X1488.0.0)
                                            ,TRUE ~ as.numeric(NA))
                ,X1488.1.0.recoded=dplyr::case_when(X1488.1.0 > 0 ~ as.numeric(X1488.1.0)
                                             ,TRUE ~ as.numeric(NA))
                ,X1488.2.0.recoded=dplyr::case_when(X1488.2.0 > 0 ~ as.numeric(X1488.2.0)
                                             ,TRUE ~ as.numeric(NA))
                ,X1498.0.0.recoded=dplyr::case_when(X1498.0.0 > 0 ~ as.numeric(X1498.0.0)
                                             ,TRUE ~ as.numeric(NA))
                ,X1498.1.0.recoded=dplyr::case_when(X1498.1.0 > 0 ~ as.numeric(X1498.1.0)
                                             ,TRUE ~ as.numeric(NA))
                ,X1498.2.0.recoded=dplyr::case_when(X1498.2.0 > 0 ~ as.numeric(X1498.2.0)
                                             ,TRUE ~ as.numeric(NA))
                ,X1508.0.0.recoded=dplyr::case_when(X1508.0.0 %in% c(-3, -1) ~ as.numeric(NA)
                                             ,X1508.0.0 == 1 ~ 0
                                             ,X1508.0.0 %in% c(2,3,4) ~ 1)
                ,X1508.1.0.recoded=dplyr::case_when(X1508.1.0 %in% c(-3, -1) ~ as.numeric(NA)
                                             ,X1508.1.0 == 1 ~ 0
                                             ,X1508.1.0 %in% c(2,3,4) ~ 1)
                ,X1508.2.0.recoded=dplyr::case_when(X1508.2.0 %in% c(-3, -1) ~ as.numeric(NA)
                                             ,X1508.2.0 == 1 ~ 0
                                             ,X1508.2.0 %in% c(2,3,4) ~ 1)
                ) %>% # Note that as.numeric(NA) is the numeric version of NA
  # Use the my_rowmeans() to take the mean from 3 recoded instances while excluding NAs
  ## Place the comma that separate new column A and B right after A. Placing it in a new line won't work
  dplyr::mutate(X1488.mean= my_rowmeans( X1488.0.0.recoded
                                        ,X1488.1.0.recoded
                                        ,X1488.2.0.recoded
                                        ,na.rm=TRUE),
                X1498.mean= my_rowmeans(X1498.0.0.recoded
                                        ,X1498.1.0.recoded
                                        ,X1498.2.0.recoded
                                        ,na.rm=TRUE),
                X1508.mean= my_rowmeans( X1508.0.0.recoded
                                         ,X1508.1.0.recoded
                                         ,X1508.2.0.recoded
                                         ,na.rm=TRUE)) %>%
  # Categorise coffee drinkers as consuming caffeine (1), drinking decaffeinated coffee (i.e. no caffeine; 0), or NA (no data) using the frequency table below 
  dplyr::mutate(consume.caffeine.thru.coffee=dplyr::case_when( X1508.mean==0 ~  0
                                               ,X1508.mean==1 ~  1
                                               ,TRUE ~ as.numeric(NA))) %>%
  # Calculate caffeine consumption per day
  dplyr::mutate(caffeine.per.day= dplyr::case_when( consume.caffeine.thru.coffee==0 ~ 40*X1488.mean
                                            ,consume.caffeine.thru.coffee==1 ~ 40*X1488.mean + 75*X1498.mean
                                            ,TRUE ~ as.numeric(NA))) # dim(ukb.pheno.1488.1498.1508) 487409     33
  
# Check frequencies
as.data.frame(table(ukb.pheno.1488.1498.1508$X1508.mean,exclude=FALSE))

# Check frequencies in X1508.mean and recoding values as consume.caffeine.thru.coffee. Note combinations in Possible are simply occurrences (i.e. 1,0,0 can be <1> 1,0,0, <2> 0,1,0, <3> 0,0,1 )

#                Var1   Freq  Possible    Recoding in consume.caffeine.thru.coffee
#----------------------------------------------------------------------------------------------
# 1                 0  71,496 (0,NA,NA)   0 (no caffeine consumption thru decaffeinated coffee)
#                             (0,0,NA)    0 (no caffeine consumption thru decaffeinated coffee)
#                             (0,0,0)     0 (no caffeine consumption thru decaffeinated coffee)
# 2 0.333333333333333    209  (1,0,0)     NA (inconsistent responses set to missing)
# 3               0.5   2496  (1,0,NA)    NA (inconsistent responses set to missing; )
# 4 0.666666666666667    314  (1,1,0)     NA (inconsistent responses set to missing)
# 5                 1 303586  (1,1,1)     1 (consume caffeine)
#                             (1,1,NA)    1 (consume caffeine)
#                             (1,NA,NA)   1 (consume caffeine)
# 6               NaN 109308  (NA,NA,NA)  NA (no data available set to missing)
#----------------------------------------------------------------------------------------------

hist(ukb.pheno.1488.1498.1508$caffeine.per.day)
summary(ukb.pheno.1488.1498.1508$caffeine.per.day)

# Export phenotype for running observation associations (this doesn't exclude cannabis users)
ExportFileSpaceSeparated(data=ukb.pheno.1488.1498.1508
                         ,missing.values.as = "NA"
                         ,output.file.path = paste0(lab.stuartma.pheno.caffeine,"daily.caffeine.consumption.thru.coffee.tea"))

#------------------------------------------------------------------------------------------------------
# Subset caffeine use data from those without cannabis initiation data. Drop y from x where y matches x
#------------------------------------------------------------------------------------------------------
# Order columns in the joined data
columns.select.ordered <- c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453_0_0_recoded","caffeine.per.day")

ukb.pheno.1488.1498.1508_20453.NA <- dplyr::inner_join(x= ukb.pheno.1488.1498.1508 # dim(ukb.pheno.1488.1498.1508) 487409 34
                                                      ,y=pheno_ukb20453_NA # dim(pheno_ukb20453_NA) 333683  3
                                                      ,by=c("FID"="FID"
                                                            ,"IID"="IID")) %>%
  dplyr::select_(.dots=columns.select.ordered) # dim(ukb.pheno.1488.1498.1508_20453.NA) 333683 12

# Export the phenotype to the BOLT-LMM folder then run GWAS
output.folder.path.UKB.caffeine.NA20453 <- paste0(lab.stuartma.gwas,"/BOLT_LMM/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/")

dir.create(output.folder.path.UKB.caffeine.NA20453)

ExportFileSpaceSeparated( data = ukb.pheno.1488.1498.1508_20453.NA
                         ,missing.values.as = "NA"
                         ,output.file.path = paste0(output.folder.path.UKB.caffeine.NA20453,"phenotype"))

#-----------------------------------------------------------------------------------------------
# Combine phenotypes 
#-----------------------------------------------------------------------------------------------
# Merge first 9 files as 1 file
list.ukb.phenotypes <- list( pheno_ukb20160
                             ,ukb.pheno.1488.1498.1508[,c("FID","IID","caffeine.per.day")])

ukb.phenotypes <- plyr::join_all(list.ukb.phenotypes,by=c("FID","IID"),type = "inner") # dim(ukb.phenotypes) 487409     16


#------------------------------------------------------------------------------------------------------
# Stratify caffeine use by smoking status (never smokers and ever smokers
#------------------------------------------------------------------------------------------------------

table(pheno_ukb20160$X20160_recode)

ukb.pheno.1488.1498.1508.20160.ever.smokers <- 

########################################################################################
#--------------3456 Number of cigarettes currently smoked daily (current smokers)------#
########################################################################################
dir_ukb3456 <- paste0(dir_ukbPheno,"/ukb3456_numCigareDaily")
file_ukb3456 <- "ukb3456.phenoUtility"
pheno_ukb3456 <- read.table(paste0(dir_ukb3456,"/",file_ukb3456),sep="",header = T) # dim(pheno_ukb3456) 487409     13

hist(pheno_ukb3456$X3456.0.0) # Use this variable for running GWAS, as X3456.0.1 and X3456.0.2 are empty
unique(pheno_ukb3456$X3456.0.1) #NULL
unique(pheno_ukb3456$X3456.0.2) #NULL

########################################################################################
#------------------------------ukb 20116 smoking Status--------------------------------#
########################################################################################
dir_ukb20116 <- paste0(dir_ukbPheno,"/","ukb20116_smokingStatus")
file_ukb20116="ukb20116.phenoUtility"
pheno_ukb20116 <- read.table(paste0(dir_ukb20116,"/",file_ukb20116),header = T,sep="",stringsAsFactors = F) # dim(pheno_ukb20116) 487409 13

# Copy information about this phenotype
## From http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20116

# Instance Description
#-----------------------------------------------------------------------------------------
# 0         Initial assessment visit (2006-2010) at which participants were recruited and consent given 501,724 participants, 501,724 items
# 1         First repeat assessment visit (2012-13) 20,339 participants, 20,339 items 
# 2         Imaging visit (2014+) 21,314 participants, 21,314 items
#-----------------------------------------------------------------------------------------

# Check phenotype value frequency
table(pheno_ukb20116$X20116.0.0,exclude = NULL)
# -3      0      1      2   <NA> 
# 1986 265431 168254  51208    530 

# Recode each of 3 instances
table(pheno_ukb20116$X20116.1.0,exclude = NULL)
# -3      0      1      2   <NA> 
# 50  12065   7159    919 467216

table(pheno_ukb20116$X20116.2.0,exclude = NULL)
# -3      0      1      2   <NA> 
# 32   7770   4469    562 474576 

# Code  Definition            NewCode everSmoker
#-------------------------------------------------
# -3	  Prefer not to answer  NA      NA
# 0	    Never                 NA      0
# 1	    Previous              0       1
# 2	    Current               1       1
#---------------------------------------------------

# Set values to NA if old values are -3 or 0. 
# Set values to 1 if old values are 2. 
# Set values to 0 if old values are 1. 

## process instance 0
pheno_ukb20116$X20116.0.0_recode <- ifelse(pheno_ukb20116$X20116.0.0 %in% c(-3,0), NA
                                           ,ifelse(pheno_ukb20116$X20116.0.0==1,0
                                                   ,1))
table(pheno_ukb20116$X20116.0.0_recode,pheno_ukb20116$X20116.0.0,useNA = "ifany")

## process instance 1
pheno_ukb20116$X20116.1.0_recode <- ifelse(pheno_ukb20116$X20116.1.0 %in% c(-3,0), NA
                                           ,ifelse(pheno_ukb20116$X20116.1.0==1,0
                                                   ,1))
table(pheno_ukb20116$X20116.1.0_recode,pheno_ukb20116$X20116.1.0,useNA = "ifany")

## process instance 2
pheno_ukb20116$X20116.2.0_recode <- ifelse(pheno_ukb20116$X20116.2.0 %in% c(-3,0), NA
                                           ,ifelse(pheno_ukb20116$X20116.2.0==1,0
                                                   ,1))
table(pheno_ukb20116$X20116.2.0_recode,pheno_ukb20116$X20116.2.0,useNA = "ifany")

# Take an average of the 3 instances
pheno_ukb20116$X20116_recode_avg= rowMeans(pheno_ukb20116[,c(14:16)],na.rm = TRUE)

as.data.frame(table(pheno_ukb20116$X20116_recode_avg,exclude=NULL))

#                Var1   Freq  PossibleA Recode Meaning  
#-------------------------------------------------------------------------------
# 1                 0 168472  (0,NA,NA)
#                             (0,0,NA)
#                             (0,0,0)   1     Never smokers
# 2 0.333333333333333    100  (1,0,0)   NA
# 3               0.5    818  (1,0,NA)  NA
# 4 0.666666666666667     48  (1,1,0)
#                             (1,1,NA)  NA
# 5                 1  50455  (1,1,1)   
#                             (1,1,NA)
#                             (1,NA,NA) 2     Ever smokers
# 6               NaN 267516  (NA,NA,NA)
#--------------------------------------------------------------------------------
pheno_ukb20116$X20116_recodeFinal <- ifelse(pheno_ukb20116$X20116_recode_avg==0,1
                                            ,ifelse(pheno_ukb20116$X20116_recode_avg==1,2
                                                    ,NA))

table(pheno_ukb20116$X20116_recodeFinal,exclude=NULL)
#     1      2   <NA> 
# 168,472  50,455 268482 

# Export processed file back to the pheno folder
write.table(pheno_ukb20116
            ,file =paste0(dir_ukb20116,"/",file_ukb20116,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

########################################################################################
#-------------------------------ukb 20160 ever smoked----------------------------------#
# Never smokers versus ever smokers by 20160 (193914, 288789) is very different from by 20116 (168472  50455)
# URL: http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20160
# Note sample sizes vary considerably from instance 0 (N=499,638), instance 1 (20,286), to instance 2 (42,724) and instance 3 (826). Determining cases or controls using an average can be difficult and misleading
########################################################################################
# Raw phenotype
file_ukb20160 <- paste0(dir_ukbPheno,"/","ukb20160_everSmoked/ukb20160.phenoUtility")
pheno_ukb20160 <- read.table(file_ukb20160,sep="",header = T) # dim(pheno_ukb20160) 487409     13

# Take an average of 3 visits
pheno_ukb20160$X20160_mean= rowMeans(pheno_ukb20160[,c(11:13)],na.rm = TRUE)
table(pheno_ukb20160$X20160_mean,exclude = NULL)

# Recode the phenotype values above according to the consistency in responses from the 3 instances. If the 3 instances contain
## 1 yes and zero No, 2 yes and zero No, or 3 yes, then recoded value set to 2
## 1 No and zero Yes, 2 No and zero Yes, or 3 No, then recoded value set to 1
## both yes and no, then recoded value set to NA (i.e. drop)

# mean  PossibleA   Count   Recode  Meaning
#-------------------------------------------------------------------------------
# 0     (0,0,0)     193914  1       Never smokers
#       (0,0,NA)
#       (0,NA,NA)
# 0.33  (1,0,0)     213     NA

# 0.5   (1,0,NA)    1905    NA
# 0.66  (1,1,0)     191     NA
# 1     (1,1,1)     288789  2       Ever smokers
#       (1,1,NA)
#       (1,NA,NA)
# NaN   (NA,8NA,NA)         NA
#--------------------------------------------------------------------------------

# Recode the value based on variable pheno_ukb20160$X20160_mean
pheno_ukb20160$X20160_recode <- pheno_ukb20160$X20160_mean

# Set consistent No (value=0) to 1 (never smokers)
pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean==0] =1
# Set consistent Yes (value=1) to 2 (ever smokers)
pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean==1] =2
# Set inconsisent responses to NA
pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean<1 & pheno_ukb20160$X20160_mean >0]= NA 
#pheno_ukb20160$X20160_recode[pheno_ukb20160$X20160_mean==NaN]=NA

table(pheno_ukb20160$X20160_recode,exclude=NULL)
#      1      2   <NA>    NaN 
# 193914 288789   2309   2397

# Export processed file back to the pheno folder
write.table(pheno_ukb20160
            ,file =paste0(file_ukb20160,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)


########################################################################################
#--------------------------ukb 3436 Age started smoking in current smokers ------------#
########################################################################################
dir_ukb3436 <- paste0(dir_ukbPheno,"/","ukb3436_ageStartedSmokingInCurrentSmokers")
file_ukb3436="ukb3436.phenoUtility"
pheno_ukb3436 <- read.table(paste0(dir_ukb3436,"/",file_ukb3436),header = TRUE,sep="",stringsAsFactors = F)

# Code                NewCode
#----------------------------------
# -1	Do not know           NA
# -3	Prefer not to answer  NA
#----------------------------------
class(pheno_ukb3436$X3436.0.0) # integer
unique(pheno_ukb3436$X3436.0.0) # 65 distinct values, consistent with http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=3436 62 distinct values if exclude NA, -3 and -1

unique(pheno_ukb3436$X3436.1.0)

unique(pheno_ukb3436$X3436.2.0)

# Set -3, -1 to NA; else kept same
## Process instance 0
pheno_ukb3436$X3436.0.0_recode <- ifelse(pheno_ukb3436$X3436.0.0 %in% c(-3,-1),NA
                                      ,pheno_ukb3436$X3436.0.0)
unique(pheno_ukb3436$X3436.0.0_recode) # -3, -1 disappear

## Process instance 1
pheno_ukb3436$X3436.1.0_recode <- ifelse(pheno_ukb3436$X3436.1.0 %in% c(-3,-1),NA
                                      ,pheno_ukb3436$X3436.1.0)
unique(pheno_ukb3436$X3436.1.0_recode) # -3, -1 disappear

## Process instance 2
pheno_ukb3436$X3436.2.0_recode <- ifelse(pheno_ukb3436$X3436.2.0 %in% c(-3,-1),NA
                                      ,pheno_ukb3436$X3436.2.0)
unique(pheno_ukb3436$X3436.2.0_recode) # -3, -1 disappear

# Take an average of the 3 recoded variables
pheno_ukb3436$X3436_recodeMean= rowMeans(pheno_ukb3436[,c(14:16)],na.rm = TRUE)

hist(pheno_ukb3436$X3436_recodeMean)

# Export processed file back to the pheno folder
write.table(pheno_ukb3436
            ,file =paste0(dir_ukb3436,"/",file_ukb3436,".recoded")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

########################################################################################
#--------------------1559 number of standard drinks per week---------------------------# 
########################################################################################
dir_ukb1559="/reference/data/UKBB_500k/versions/lab_stuartma/pheno/alcohol"
file_ukb1559="alcohol.recoded.weeklyunits.full.pheno"
pheno_ukb1559 <- read.table(paste0(dir_ukb1559,"/",file_ukb1559),header = TRUE,sep="",stringsAsFactors = F) #dim(pheno_ukb1559) 487409     13

# Rename the wanted variable as it is too long for HPC-Utility
pheno_ukb1559$NSDPW <- pheno_ukb1559$alcrecoded_weeklyunits

hist(pheno_ukb1559$NSDPW)
unique(pheno_ukb1559$NSDPW)

# Use variable pheno_ukb1559$NSDPW as the phenotype for running GWAS. No processing for this variable as it is calcualted by JueSheng Ong

# Copy the analysed file to my own phenotype folder
#output="/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb1559_numberStDrinksPerWeek"
output <- paste0(dir_ukbPheno,"/ukb1559_numberStDrinksPerWeek")
write.table(pheno_ukb1559
            ,file =paste0(output,"/",file_ukb1559)
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

########################################################################################
#--------------3456 Number of cigarettes currently smoked daily (current smokers)------# 
#-------------- excluding those with data in 20453 ever using cannabis ----------------#
########################################################################################

# Left join pheno_ukb20453_NA (left table) and pheno_ukb3456 (right table)
pheno_ukb3456_NA20453 <-  merge(x=pheno_ukb20453_NA, y=pheno_ukb3456, by="IID", all.x = TRUE) # 333683 obs. of  15 variables

# Check the 3 instances of X3456
table(pheno_ukb3456_NA20453$X3456.0.0) # values ranged from -10 to 140
table(pheno_ukb3456_NA20453$X3456.1.0) # values ranged from -3 to 60
table(pheno_ukb3456_NA20453$X3456.2.0) # values ranged from -1 to 30

# Take an average of X3456.0.0, X3456.1.0, X3456.2.0
## First set negative values to NA; other values as they are originally
pheno_ukb3456_NA20453$X3456.0.0_recode <-  ifelse(pheno_ukb3456_NA20453$X3456.0.0 < 0, NA
                                               ,pheno_ukb3456_NA20453$X3456.0.0)

pheno_ukb3456_NA20453$X3456.1.0_recode <-  ifelse(pheno_ukb3456_NA20453$X3456.1.0 < 0, NA
                                               ,pheno_ukb3456_NA20453$X3456.1.0)

pheno_ukb3456_NA20453$X3456.2.0_recode <-  ifelse(pheno_ukb3456_NA20453$X3456.2.0 < 0, NA
                                               ,pheno_ukb3456_NA20453$X3456.2.0)  
# Take an average of 3 visits
pheno_ukb3456_NA20453$X3456_mean <-  rowMeans(pheno_ukb3456_NA20453[,c(16:18)],na.rm = TRUE)

# Rename columns
pheno_ukb3456_NA20453$FID <- pheno_ukb3456_NA20453$FID.x

# Reorder columns
columns_want_in_this_order <- c("FID","IID","missing","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X3456.0.0","X3456.1.0","X3456.2.0","X3456.0.0_recode","X3456.1.0_recode","X3456.2.0_recode","X3456_mean")

pheno_ukb3456_NA20453_ordered <- pheno_ukb3456_NA20453 %>% select_(.dots=columns_want_in_this_order)

# Export the phenotype data to the 3456 folder as a backup
write.table(pheno_ukb3456_NA20453_ordered
            ,file=paste0(dir_ukb3456,"/ukb3456_IID_NA_in_20453")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputFolderPath_ukb3456_NA20453 <- paste0(lab.stuartma.gwas,"/BOLT_LMM/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/phenotype")

write.table(pheno_ukb3456_NA20453_ordered
            ,file=paste0(outputFolderPath_ukb3456_NA20453,"/ukb3456_IID_NA_in_20453")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

########################################################################################
#--------------UKB complete_alcohol_unitsweekly
#-------------- excluding those with data in 20453 ever using cannabis ----------------#
########################################################################################
# This computed variable combined weekly consumption in heavy drinkers (ukb1558=4,5,6) and monthly consumption in light drinkers (ukb1558=1,2,3)

# Location of the phenotype
fileName_ukb_alcoholUnitsWeekely <- "total_alcohol_unitsweekly.combined.pheno"

# Import the weekly alcohol consumption file
file_alcoholUnitsWeekely <- read.table(paste0(lab_stuartma_pheno_alcohol,fileName_ukb_alcoholUnitsWeekely),sep=" ",header = T) # 487409 obs. of  32 variables

# Exclude IIDs with data for 20453 ever using cannabis
## Left join pheno_ukb20453_NA (left table) and file_alcoholUnitsWeekely (right table)
pheno_ukbAlcUnitsWeekly_NA20453 <-  merge(x=pheno_ukb20453_NA
                                       ,y=file_alcoholUnitsWeekely
                                       , by=c("FID","IID")
                                       , all.x = TRUE) # 333683 obs. of  33 variables

# Reorder columns
columns_ordered <- c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453.0.0","complete_alcohol_unitsweekly")
pheno_ukbAlcUnitsWeekly_NA20453_ordered <- pheno_ukbAlcUnitsWeekly_NA20453 %>% 
                                            select_(.dots=columns_ordered)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputFolderPath_ESDPW_NA20453 <- paste0(lab.stuartma.gwas,"/BOLT_LMM/UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis")
dir.create(outputFolderPath_ESDPW_NA20453)

write.table(pheno_ukbAlcUnitsWeekly_NA20453_ordered
            ,file=paste0(outputFolderPath_ESDPW_NA20453,"/phenotype")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

########################################################################################
#--------------UKB 20161 Pack years of smoking ------------------------------------#
#-------------- excluding those with data in 20453 ever using cannabis ------------#
# ULR: http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20161
########################################################################################
# Raw phenotype location
file_ukb20161 <- paste0(lab_stuartma_pheno_smoking,"pack_years.pheno")
pheno_ukb20161 <- read.table(file_ukb20161,sep=" ",header = T,stringsAsFactors = F)
# JueSheng has merged the 3 instances to the column above. Values from non smokers set to 0 so there are not too many NAs

# Exclude IIDs with data for 20453 ever using cannabis
## Left join pheno_ukb20453_NA (left table) and pheno_ukb20161 (right table)
pheno_ukb20161_NA20453 <-  merge(x=pheno_ukb20453_NA 
                              ,y=pheno_ukb20161
                              , by=c("FID","IID")
                              , all.x = TRUE) # 333683 obs. of  17 variables

# Reorder columns
columns_ordered <- c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453.0.0","merged_pack_years_20161")

pheno_ukb20161_NA20453_ordered <- pheno_ukb20161_NA20453 %>% select_(.dots=columns_ordered)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputFolderPath_ukb20161_NA20453 <- paste0(lab.stuartma.gwas,"/BOLT_LMM/UKB20161-pack-years-of-smoking_IID-NA-in-UKB20453-everUsedCannabis")

dir.create(outputFolderPath_ukb20161_NA20453)

write.table(pheno_ukb20161_NA20453_ordered
            ,file=paste0(outputFolderPath_ukb20161_NA20453,"/phenotype")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

########################################################################################
#--------------UKB cups of coffee per day (Data-Field 1498 coffee intak?)--------------#
#-------------- excluding those with data in 20453 ever using cannabis ----------------#
########################################################################################
# Raw phenotype location
file_ukb_CCPD <- paste0(lab_stuartma_pheno,"coffee/new_coffee.pheno")
pheno_ukb_CCPD <- read.table(file_ukb_CCPD,sep=" ",header = T,stringsAsFactors = F)

# Exclude IIDs with data for 20453 ever using cannabis
## Left join pheno_ukb20453_NA (left table) and pheno_ukb_CCPD (right table)
pheno_ukb_CCPD_NA20453 <-  merge(x=pheno_ukb20453_NA
                              ,y=pheno_ukb_CCPD
                              , by=c("FID","IID")
                              , all.x = TRUE) # 333683 obs. of  21 variables

# Reorder columns
columns_ordered <- c("FID","IID","missing","batch","kinship","exclude_kinship","excess_relative","age","sex","white.British","X20453.0.0","all_coffee_cpd")

pheno_ukb_CCPD_NA20453_ordered <- pheno_ukb_CCPD_NA20453 %>% select_(.dots=columns_ordered)

# Export the phenotype data to the BOLT-LMM for running GWAS using HPC_Utility.jar
outputFolderPath_ukb_CCPD_NA20453 <- paste0(lab.stuartma.gwas,"/BOLT_LMM/UKB-cups-of-coffee-per-day_IID-NA-in-UKB20453-everUsedCannabis")

dir.create(outputFolderPath_ukb_CCPD_NA20453)

write.table(pheno_ukb_CCPD_NA20453_ordered
            ,file=paste0(outputFolderPath_ukb_CCPD_NA20453,"/phenotype")
            ,sep=" "
            ,na="NA"
            ,col.names = TRUE
            ,row.names = FALSE
            ,quote=FALSE)

# file.copy("/mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step00-00_recode_phenotype_for_GWAS.R","/mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step07-01_extract-UKB-phenotypes.R")

#---------------------------------------------------------------------------------------#
#--------------------------- This is the end of this program----------------------------#
#---------------------------------------------------------------------------------------#

