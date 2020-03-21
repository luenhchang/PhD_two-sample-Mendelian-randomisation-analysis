##################################################################################
# Filename: MR_step07_observational-analysis_effect-exposure-trait-on-outcome-trait_UKB-phenotype-data.R
# Programmer: Chang
# Purpose: (1) Examine observational assocation between exposure traits and outcome traits using UKB phenotype data (2) calculate %females, age from participants with non-missing data in outcome of interest
# Date created: 20180819
# Dependent: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step00-01_recode-UKB-phenotypes-for-running-GWAS.R
# Note: (1) UKB GWAS conducted at prior steps use phenotypes from participants who didn't respond to ever taken cannabis. Don't use those phenotypes here, (2) binary outcome variables should be 1 for cases, 0 for control
#-----------------------------------------------------------------------------------------

# Type 	File
#------------------------------------------------------------------------------------------------
# Input	paste0(stuartma_pheno,"coffee/additional_coffee_lifestyle.pheno")
# Input	paste0(stuartmaDir,"collab/alcohol_BCAC/total_alcohol_unitsweekly.combined.pheno")
# Input	paste0(stuartma_pheno,"coffee/new_coffee.pheno")
# Input paste0(stuartma.pheno.caffeine,"daily.caffeine.consumption.thru.coffee.tea")
# Input	paste0(stuartma_pheno,"/may18_baseline_extra_covariate")
# Input	paste0(dir_ukbPheno,"/ukb20116_smokingStatus/ukb20116.phenoUtility.recoded")
# Input	paste0(dir_ukbPheno,"/ukb20160_everSmoked/ukb20160.phenoUtility.recoded")
# Input	paste0(stuartma_pheno,"smoking/pack_years.pheno")
# Input	paste0(dir_ukbPheno,"/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded")
# Input	paste0(dir_ukbPheno,"/ukb3436_ageStartedSmokingInCurrentSmokers/ukb3436.phenoUtility.recoded")
# Input	paste0(dir_ukbPheno,"/ukb3456_numCigareDaily/ukb3456.phenoUtility")
# Input	paste0(stuartmaDir,"exclude_related/notWhite.id")
# Input	paste0(stuartmaDir,"relatives_working/alltrait_rela")
# Input	paste0(loc.obs.assoc,"exposure-outcome-pairs-used-in-MR-analyses.tsv")

# Outpu paste0(loc.obs.assoc,"results_binary-logistic-regression_linear-regression_full-parameters.tsv") 
# Outpu paste0(loc.obs.assoc,"sample-sizes_phenotypic-associations.tsv")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190930  Included PC9, PC10 in the regression. Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv
# 20190928  Exported paste0(loc.obs.assoc,"sample-sizes_phenotypic-associations.tsv")
# 20190904  Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv
# 20190903  Added ever smoked as a predictor. Standardised pack years of smoking but haven't used it in the regression. Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv
# 20190815  Copied number of observations analysed 13917 to manuscript. This number is the same from iteration 1 to 4. So no output file
# 20190813  Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv
# 20190807  Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv
# 20190806  Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv
# 20190805  Exported results_binary-logistic-regression_linear-regression_full-parameters.tsv (replace old file results_binary-logistic-regression_linear-regression.tsv)
# 20190713  Exported results_binary-logistic-regression_linear-regression.tsv
# 20190713  Rerun regression analyses adding PC1:PC8 into covariates. all modelling results combined and exported.
# 20190705  Run logistic regression of coffee, alcohol... on cannabis with 2 SNPs as covariates. Run linear regression of cannabis, alcohol... on coffee with 2 SNPs as covariates
# 20190610  Calculated (1) number of participants with missing data in their cannabis initiation 20453, (2) N participants with missing 20453 and non-missing CCPD, (3) N participants with missing 20453 and non-missing CPD, (4) N participants with missing 20453 and non-missing ESDPW, (5) N participants with missing 20453 and non-missing PYOS, (6) % females in (1) to (5)

# 20190418 Exported the 1 file above
# 20190415 Exported the 1 file above
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Folder locations under Stuart MacGregor's lab
#-----------------------------------------------------------------------------------------
stuartmaDir <- "/reference/data/UKBB_500k/versions/lab_stuartma/"
stuartma_pheno <- paste0(stuartmaDir,"pheno/")
stuartma.pheno.coffee <- paste0(stuartma_pheno,"coffee/")
stuartma.pheno.caffeine <- paste0(stuartma.pheno.coffee,"ukb_1488-teaIntake_1498-coffeeIntake_1508-coffeeType/")

#-----------------------------------------------------------------------------------------
# Folder locations under my home directory
#-----------------------------------------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
dir_ukbPheno <- paste0(homeDir,"data/UKBiobank_phenotype/")

#-----------------------------------------------------------------------------------------
# Folder locations under my working directory
#-----------------------------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/") # location of outcome data
loc.obs.assoc <- paste0(locMR,"observational-associations/")
#dir.create(loc.obs.assoc)

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#-------------------------------------------------------------------------------------#
# Import a list of exposures and outcomes used in MR analyses
# Match predictors and outcomes from the UKB to the MR exposures and outcomes
#-------------------------------------------------------------------------------------#
ImportATabSeparatedFile(input.file.path = paste0(loc.obs.assoc,"exposure-outcome-pairs-used-in-MR-analyses.tsv")
                        ,data.name = "MR.exposure.outcome.traits") # dim(MR.exposure.outcome.traits) 23 4

#-------------------------------------------------------------------------------------#
#------File paths of phenotypes, covariates, non-white, relatives---------------------#
#--- Copy recoded values and definition under each file path
#-------------------------------------------------------------------------------------#
# Education attainment
file.path.ukb.pheno.edu.att <- paste0(stuartma_pheno,"coffee/additional_coffee_lifestyle.pheno")

# Estimated standard drinks per week
file.path.ukb.pheno.ESDPW <- paste0(stuartmaDir,"collab/alcohol_BCAC/total_alcohol_unitsweekly.combined.pheno")

# Caffeine consumed per day (mg)
file.path.ukb.pheno.CCPD <- paste0(stuartma.pheno.caffeine,"daily.caffeine.consumption.thru.coffee.tea")

# Covariates
file.path.ukb.pheno.covariates <- paste0(stuartma_pheno,"may18_baseline_extra_covariate")

# Ever smoked
file.path.ukb.pheno.20160 <- paste0(dir_ukbPheno,"ukb20160_everSmoked/ukb20160.phenoUtility.recoded")

# Pack years of smoking
file.path.ukb.pheno.20161 <- paste0(stuartma_pheno,"smoking/pack_years.pheno")

# Ever taken cannabis
file.path.ukb.pheno.20453 <- paste0(dir_ukbPheno,"/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded")

# Age at starting smoking in current smokers
file.path.ukb.pheno.3436 <- paste0(dir_ukbPheno,"/ukb3436_ageStartedSmokingInCurrentSmokers/ukb3436.phenoUtility.recoded")

# Cigarettes per day
file.path.ukb.pheno.3456 <- paste0(dir_ukbPheno,"/ukb3456_numCigareDaily/ukb3456.phenoUtility")

# IDs of non-caucasian
file.path.ukb.non.white <- paste0(stuartmaDir,"exclude_related/notWhite.id")

# Relatives of UKB participants
file.path.ukb.relatives <- paste0(stuartmaDir,"relatives_working/alltrait_rela")

#-------------------------------------------------------------------------------------#
#-----------------Import UKB phenotype files -----------------------------------------#
#-------------------------------------------------------------------------------------#
common.columns <- c("FID","IID")

# UKB data-field 6138 educational attainment. 2 new variables derived from qualification_edu_6138. only 1 is in use
# Coding	Meaning
#------------------------------------------------------------------
# 1	      College or University degree
# 2	      Advanced Level(A levels)/AS levels or equivalent (AS and A levels are at level 3 on the RQF)
# 3	      O levels/General Certificate of Secondary Education(GCSEs) or equivalent (GCSEs are at levels 1 and 2 on the RQF)
# 4	      General Certificate of Education (CSEs) or equivalent
# 5	      National Vocational Qualifications (NVQ) or Higher National Diplomas (HND) or Higher National Certificates (HNC) or equivalent (HNCs and HNDs are at level 5 on the NQF)
# 6	      Other professional qualifications eg: nursing, teaching
#-----------------------------------------------------------------

# Coding	Recode
#------------------------------------------------------------------
# 1	      5
# 2	      4
# 3	      3
# 4	      2
# 5	      1
# 6	      NA
#-----------------------------------------------------------------

ukb.pheno.edu.att <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.edu.att
                                               ,data.name = "ukb.pheno.edu.att") %>%
  dplyr::select_(.dots=c(common.columns,"qualification_edu_6138")) # dim(ukb.pheno.edu.att) 487409      3

# Frequencies of original values
table(ukb.pheno.edu.att$qualification_edu_6138,exclude = FALSE)
#      1      2      3      4      5      6   <NA> 
# 158739  54515 103159  26102  31779  24550  88565

# Collapse original coding into a binary variable (value=1 or value=2 to 6)
ukb.pheno.edu.att$EA.6138.uniCollegeDeg <- ifelse(ukb.pheno.edu.att$qualification_edu_6138==1, 1,
                                                  ifelse(ukb.pheno.edu.att$qualification_edu_6138 %in% c(2:6),0, NA)
                                                  )

table(ukb.pheno.edu.att$EA.6138.uniCollegeDeg, exclude = FALSE)
#      0      1   <NA> 
# 240105 158739  88565

# Recode qualification_edu_6138 
ukb.pheno.edu.att <- ukb.pheno.edu.att %>% dplyr::mutate(quali.edu.6138.recoded=dplyr::recode(qualification_edu_6138
                                                                                              ,`1`=5
                                                                                              ,`2`=4
                                                                                              ,`3`=3
                                                                                              ,`4`=2
                                                                                              ,`5`=1))
table(ukb.pheno.edu.att$quali.edu.6138.recoded, exclude = F)
#     1      2      3      4      5   <NA> 
# 31779  26102 103159  54515 158739 113115

# Standardise the recoded variable to generate Z scores
ukb.pheno.edu.att$quali.edu.6138.recoded.z.score <- scale(ukb.pheno.edu.att$quali.edu.6138.recoded
                                                          , center = TRUE
                                                          , scale = TRUE)

ukb.pheno.ESDPW <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.ESDPW
                                             ,data.name = "ukb.pheno.ESDPW") %>%
  dplyr::select_(.dots=c(common.columns,"complete_alcohol_unitsweekly")) # dim(ukb.pheno.ESDPW) 487409 3

ukb.pheno.CCPD <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.CCPD
                                            ,data.name = "ukb.pheno.CCPD") %>%
  dplyr::select_(.dots=c(common.columns,"caffeine.per.day")) # dim(ukb.pheno.CCPD) 487409      3

# Covariates (fit these as covariates as minimum: age, inferred.sex, TDI [townsend deprivation index, a measure for poverty], maybe BMI (probably not, as it is heritable), overall_health_rating)
ukb.pheno.covar <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.covariates
                                             ,data.name = "ukb.pheno.covar") %>%
  dplyr::select_(.dots=c(common.columns,"inferred.sex","age",paste0("PC",c(1:10)),"TDI","overall_health_rating","rs2472297", "rs4410790")) # dim(ukb.pheno.covar) 487409      18

# Smoking status (for matching GSCAN smoking cessation)
ukb.pheno.20116 <- ImportASpaceSeparatedFile(input.file.path = "/mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20116_smokingStatus/ukb20116.phenoUtility.recoded"
                                             ,data.name = "ukb.pheno.20116") %>% 
  dplyr::select_(.dots=c(common.columns,"X20116_recodeFinal")) %>%
  dplyr::mutate(X20116_recodeFinal_0_1=dplyr::case_when( X20116_recodeFinal==1 ~ 0
                                                        ,X20116_recodeFinal==2 ~ 1
                                                        , TRUE ~ as.numeric(NA)))#

# Previous smokers coded as 1, current smokers coded as 2 in project PRS_UKB_201711
table(ukb.pheno.20116$X20116_recodeFinal, exclude = FALSE)
#      1      2   <NA> 
# 168472  50455 268482

# For logistic regression, outcome variables are coded as 0, 1
table(ukb.pheno.20116$X20116_recodeFinal_0_1, exclude = FALSE)
#      0      1   <NA> 
# 168472  50455 268482 

# Ever smoked (ever smokers versus never smokers)
ukb.pheno.20160 <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.20160
                                             ,data.name = "ukb.pheno.20160") %>%
  dplyr::select_(.dots=c(common.columns,"X20160_recode")) %>% 
  # Recode 2 (ever smoked) as 1, 1 (never smoked) as 0
  dplyr::mutate(X20160_recode=dplyr::recode(X20160_recode
                                     ,`2`=1
                                     ,`1`=0)) # dim(ukb.pheno.20160) 487409      3
table(ukb.pheno.20160$X20160_recode, exclude = FALSE)

# Pack years of smoking.
## Why is this variable set to zero in never smokers (X20160_recode==0)? It is supposed to be NA.
ukb.pheno.20161 <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.20161
                                             ,data.name = "ukb.pheno.20161") %>%
  dplyr::select_(.dots=c(common.columns,"merged_pack_years_20161")) # dim(ukb.pheno.20161) 487409      3

# Calculate mean and SD of pack years of smoking
mean(ukb.pheno.20161$merged_pack_years_20161, na.rm = TRUE) # 8.139645
sd(ukb.pheno.20161$merged_pack_years_20161, na.rm = TRUE) # 15.58431

# Standardise pack years of smoking
ukb.pheno.20161$merged.pack.years.20161.z.scores <- scale(ukb.pheno.20161$merged_pack_years_20161
                                                          ,center = TRUE
                                                          ,scale = TRUE)

ukb.pheno.20453 <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.20453
                                             ,data.name = "ukb.pheno.20453") %>%
  dplyr::select_(.dots=c(common.columns,"X20453_0_0_recoded")) # dim(ukb.pheno.20453) 487184      3

ukb.pheno.3436 <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.3436
                                            ,data.name = "ukb.pheno.3436") %>%
  dplyr::select_(.dots=c(common.columns,"X3436_recodeMean")) # dim(ukb.pheno.3436) 487409      3

ukb.pheno.3456 <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.pheno.3456
                                            ,data.name = "ukb.pheno.3456") %>%
  dplyr::select_(.dots=c(common.columns,"X3456.0.0")) # dim(ukb.pheno.3456) 487409      3
                         
ukb.non.white <-   ImportASpaceSeparatedFile(input.file.path = file.path.ukb.non.white
                                             ,data.name = "ukb.non.white") # dim(ukb.non.white) 48538     1
colnames(ukb.non.white) <- "IID"

ukb.relatives <- ImportASpaceSeparatedFile(input.file.path = file.path.ukb.relatives
                                           ,data.name = "ukb.relatives") %>%
  dplyr::select_(.dots=c("ID1","ID2")) # dim(ukb.relatives) 107160      2

#-----------------------------------------------------------------------------------------------
# Count number of people with non-missing data in cannabis initiation 20453
#-----------------------------------------------------------------------------------------------
table(ukb.pheno.20453$X20453_0_0_recoded,exclude = FALSE)
#      0      1   <NA> 
# 119694  33807 333908
119694+33807 # 153501

#-----------------------------------------------------------------------------------------------
# Combine phenotype and covariate files. Exclude participants based on non.white and relative files
#-----------------------------------------------------------------------------------------------
# Merge first 9 files as 1 file
list.ukb.phenos.covar <- list( ukb.pheno.edu.att
                               ,ukb.pheno.ESDPW
                               ,ukb.pheno.CCPD
                               ,ukb.pheno.covar
                               ,ukb.pheno.20116
                               ,ukb.pheno.20160 #ukb.pheno.20116
                               ,ukb.pheno.20161
                               ,ukb.pheno.20453
                               ,ukb.pheno.3436
                               ,ukb.pheno.3456)

ukb.phenos.covar <- plyr::join_all(list.ukb.phenos.covar,by=c("FID","IID"),type = "inner") # dim(ukb.phenos.covar) 487409     26

# Exclude non-white people 
ukb.phenos.covar.nonwhite.removed <- ukb.phenos.covar %>% dplyr::anti_join(ukb.non.white) # dim(ukb.phenos.covar.nonwhite.removed) 438679     27

# Exclude related individuals from the participants based on ID1 or ID2
ukb.phenos.covar.nonwhite.related.ID1.removed <- dplyr::anti_join(ukb.phenos.covar.nonwhite.removed
                                                                  ,ukb.relatives
                                                                  ,by=c("IID" = "ID1")) # dim(ukb.phenos.covar.nonwhite.related.ID1.removed) 361223     30

ukb.phenos.covar.nonwhite.related.ID2.removed <- dplyr::anti_join(ukb.phenos.covar.nonwhite.removed
                                                                  ,ukb.relatives
                                                                  ,by=c("IID" = "ID2")) # dim(ukb.phenos.covar.nonwhite.related.ID2.removed) 361283     32

#-------------------------------------------------------------------------------------
# age and sex
## inferred.sex=0 for males, inferred.sex=1 for females 
#-------------------------------------------------------------------------------------
# Demongraphic of participants used
d <- ukb.phenos.covar.nonwhite.related.ID2.removed

# Age mean, standard devidation
round(mean(d$age, na.rm = TRUE), digits = 2) # 56.79
round(sd(d$age, na.rm = TRUE), digits = 2) # 8

#-------------------------------------------------------------------------------------
# Model 1 to 5: Run a mutiple logistic regression of binary outcome or mutiple linear regression of continuous outcome
## Outcomes: CI, SI, CCPD, PYOS and ESDPW. Excluded ASSCS and CPD because of errors in every regression with these two as outcome or predictor variables
#-------------------------------------------------------------------------------------
# Set up parts of equations
equation.covariates <- "age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8+ PC9+ PC10 + TDI + overall_health_rating + factor(inferred.sex) + quali.edu.6138.recoded.z.score" # 

# Group variables that are used as outcome or predictor
obs.predictor.outcome.matched.MR <- MR.exposure.outcome.traits %>% 
  dplyr::mutate(outcome.var=dplyr::case_when(  outcome.trait=="AI" ~ "X3436_recodeMean"
                                               ,outcome.trait=="CPD" ~ "X3456.0.0"
                                               ,outcome.trait=="SC" ~ "X20116_recodeFinal_0_1"
                                               ,outcome.trait=="SI" ~ "X20160_recode"
                                               ,outcome.trait=="CI" ~ "X20453_0_0_recoded"
                                               ,outcome.trait=="caffeine" ~ "caffeine.per.day"
                                               ,outcome.trait=="ESDPW" ~ "complete_alcohol_unitsweekly"
                                               ,outcome.trait=="PYOS" ~ "merged_pack_years_20161")
                ,outcome.type=dplyr::case_when(outcome.trait %in% c("SI","CI") ~ "binary"
                                               ,TRUE~ "continuous")
                ,predictor.var=dplyr::case_when(exposure.trait=="SI" ~ "X20160_recode"
                                                ,exposure.trait=="CI" ~ "X20453_0_0_recoded"
                                                ,exposure.trait=="caffeine" ~ "caffeine.per.day"
                                                ,exposure.trait=="ESDPW" ~ "complete_alcohol_unitsweekly")
                ,predictor.type=dplyr::case_when(exposure.trait %in% c("SI","CI") ~ "binary"
                                                 ,TRUE~ "continuous")) %>% 
  # Remove unreasonable pairs
  dplyr::filter(!(outcome.trait=="CPD" & exposure.trait=="SI")) %>%
  dplyr::rename(predictor.label=exposure.trait, outcome.label=outcome.trait) %>%
  dplyr::select(outcome.var,outcome.label,outcome.type,predictor.var,predictor.label,predictor.type) %>%
  dplyr::arrange(outcome.label, predictor.label) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(formula= paste0(outcome.var, "~", predictor.var, "+", equation.covariates)) # dim(obs.predictor.outcome.matched.MR) 21 7

# Create an empty data.frame for appending results of all iterations
base.logis.coeff <- data.frame()
base.logis.sample.size <- data.frame()
base.linea.coeff <- data.frame()
base.linea.sample.size <- data.frame()

# Loop through each row 
count <- 0

data <- ukb.phenos.covar.nonwhite.related.ID2.removed # dim(data) 361283     27

numb.obs.input.data <- nrow(data)

# Number of iterations: 22
for (i in 1:nrow(obs.predictor.outcome.matched.MR)){
  count <- count+1
  varName <- obs.predictor.outcome.matched.MR[i,"outcome.var"] # 1 used as an outcome variable, the other 5 used as predictors
  var.label <- obs.predictor.outcome.matched.MR[i,"outcome.label"]
  var.type <- obs.predictor.outcome.matched.MR[i,"outcome.type"]
  formula <- obs.predictor.outcome.matched.MR[i,"formula"]
  print(paste0("======================== iteration ",count,"============================="))
  
  tryCatch({
    if (var.type=="binary"){
      print(paste0("Running a binary logistic regression of ", varName))
      model.fit.logistic <- glm(formula= formula,data=data,family=binomial(link = "logit"))
      mod.logi.summ <- summary(model.fit.logistic) # str(mod.logi.summ) a list of 18
      
      # Get coefficients
      mod.logi.summ.coeff <- mod.logi.summ[["coefficients"]] # dim(mod.logi.summ.coeff) 12 4
      
      # Convert a matrix to a data.frame keeping its dimnames
      mod.logi.summ.coeff.df <- plyr::adply(mod.logi.summ.coeff,1,c)
      mod.logi.summ.coeff.df$X1 <- as.character(mod.logi.summ.coeff.df$X1) # dim(mod.logi.summ.coeff.df)
      
      # Rename columns
      colnames(mod.logi.summ.coeff.df) <- c("predictor","estimate","SE","z.value","p.value")
      
      # Add name of the dependent variable
      mod.logi.summ.coeff.df$dep.var <- varName
      mod.logi.summ.coeff.df$dep.var.label <- var.label

      ## Add iteration number 
      mod.logi.summ.coeff.df$iteration <- count
      
      # Append current iteration data.frame to the base data set
      base.logis.coeff <- rbind(base.logis.coeff,mod.logi.summ.coeff.df)
      
      # Calculate number of observations used in the model fitting
      numb.obs.deleted <- length(mod.logi.summ[["na.action"]])
      numb.obs.analysed <- numb.obs.input.data-numb.obs.deleted
      logis.sample.sizes <- data.frame( dep.var=varName
                                        ,dep.var.label=var.label
                                        ,iteration=count
                                        ,numb.obs.input.data=nrow(data)
                                        ,numb.obs.deleted= numb.obs.deleted
                                        ,numb.obs.analysed= numb.obs.analysed
                                        ,stringsAsFactors = F)
      # Append number of observation information to the base data
      base.logis.sample.size <- rbind(base.logis.sample.size,logis.sample.sizes)
      
      
    } else {
      print(paste0("Running a linear regression of ", varName))
      model.fit.linear <- lm(formula = formula, data=data)
      mod.line.summ <- summary(model.fit.linear)
      mod.line.summ.coeff <- mod.line.summ[["coefficients"]] # dim(mod.line.summ.coeff) 12 4 class(mod.line.summ.coeff) "matrix"
      # Convert a matrix to a data.frame keeping its dimnames
      mod.line.summ.coeff.df <- plyr::adply(mod.line.summ.coeff, 1, c) # class(mod.line.summ.coeff.df) "data.frame"
      mod.line.summ.coeff.df$X1 <- as.character(mod.line.summ.coeff.df$X1)
      # Rename columns
      colnames(mod.line.summ.coeff.df) <- c("predictor","estimate","SE","t.value","p.value")
      
      # Add name of the dependent variable
      mod.line.summ.coeff.df$dep.var <- varName
      mod.line.summ.coeff.df$dep.var.label <- var.label
      ## Add iteration number 
      mod.line.summ.coeff.df$iteration <- count
      
      # Append current iteration data.frame to the base data set
      base.linea.coeff <- rbind(base.linea.coeff, mod.line.summ.coeff.df)
      
      # Calculate number of observations used in the model fitting
      numb.obs.deleted <- length(mod.line.summ[["na.action"]])
      numb.obs.analysed <- numb.obs.input.data-numb.obs.deleted
      linea.sample.sizes <- data.frame( dep.var=varName
                                       ,dep.var.label=var.label
                                       ,iteration=count
                                       ,numb.obs.input.data=nrow(data)
                                       ,numb.obs.deleted=numb.obs.deleted
                                       ,numb.obs.analysed=numb.obs.analysed
                                       ,stringsAsFactors = F)
      # Append number of observation information to the base data
      base.linea.sample.size <- rbind(base.linea.sample.size,linea.sample.sizes)
      
    }
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

dim(base.logis.coeff) # 102 8
dim(base.logis.sample.size) # 6 6
dim(base.linea.coeff) # 255   8
dim(base.linea.sample.size) # 15  6

# Rename columns for merging purposes
base.logis.coeff2 <- base.logis.coeff %>% 
  dplyr::mutate(test.stat.name="z.value",model="logistic regression") %>%
  dplyr::rename(test.stat.value=z.value)

base.linea.coeff2 <- base.linea.coeff %>%
  dplyr::mutate(test.stat.name="t.value", model="linear regression") %>%
  dplyr::rename(test.stat.value=t.value)
  
# Merge base.logis.coeff and base.linea.coeff
coeff <- dplyr::bind_rows(base.logis.coeff2, base.linea.coeff2) # dim(coeff) 357 10

# Delete string "factor(" and ")" from column predictor
coeff$predictor <- gsub(coeff$predictor, pattern = "\\(|\\)|factor\\(|\\)1",replacement="")

# Get a list of unique outcomes
outcome.unique <- obs.predictor.outcome.matched.MR[,c("outcome.label","outcome.type")] %>% 
  dplyr::distinct() # dim(outcome.unique) 8 2

# Add labels to predictors and covariates

# Create variable type and label for covariates
covar.name.type <- data.frame(predictor.var=c("age"
                                        ,paste0("PC",c(1:10))
                                        ,"TDI"
                                        ,"overall_health_rating"
                                        ,"quali.edu.6138.recoded.z.score"
                                        ,"inferred.sex")
                              ,predictor.label=c("Age"
                                           ,paste0("PC",c(1:10))
                                           ,"TDI"
                                           ,"Overall health rating"
                                           ,"Educational attainment"
                                           ,"Sex")
                              ,predictor.type=c(rep("continuous",times=14),"binary")
                              ,stringsAsFactors = F) # dim(covar.name.type) 15 3

# Stack predictor and covariate names and labels
all.variables <- obs.predictor.outcome.matched.MR[,c( "predictor.var"
                                                      ,"predictor.label"
                                                      ,"predictor.type")] %>% 
  # Rid of duplicate rows
  dplyr::distinct() %>%
  dplyr::bind_rows(.,covar.name.type) # dim(all.variables) 19 3

caffeine.per.cup <- 0.5*(75+40)

# Add variable types and labels to predictors and covariates
coeff2 <- dplyr::left_join( coeff
                           ,all.variables
                           ,by=c("predictor"="predictor.var")) %>%
  ## Add variable types (outcome.type) to dependent variables in the coefficient data
  dplyr::left_join(.
                   ,outcome.unique
                   ,by=c("dep.var.label"="outcome.label")) %>%
  ## Rename the added columns
  dplyr::rename(dep.var.type=outcome.type, b=estimate) %>%
  ## Convert effect estimates 
  dplyr::mutate(effect.size.estimate=dplyr::case_when(
     predictor.type=="continuous" & predictor.label!="ECCPD" & dep.var.type=="continuous" ~ b
    ,predictor.type=="continuous" & predictor.label=="ECCPD" & dep.var.type=="continuous" ~ b*caffeine.per.cup
    ,predictor.type=="continuous" & predictor.label!="ECCPD" & dep.var.type=="binary" ~ exp(b)
    ,predictor.type=="continuous" & predictor.label=="ECCPD" & dep.var.type=="binary" ~ exp(b*caffeine.per.cup)
    ,predictor.type=="binary" & dep.var.type=="continuous" ~ log(2)*b
    ,predictor.type=="binary" & dep.var.type=="binary" ~ exp(log(2)*b)
                                              )
    ,effect.size.LBound=dplyr::case_when(
       predictor.type=="continuous" & predictor.label!="ECCPD" & dep.var.type=="continuous" ~ b-1.96*SE
      ,predictor.type=="continuous" & predictor.label=="ECCPD" & dep.var.type=="continuous" ~ caffeine.per.cup*(b-1.96*SE)
      ,predictor.type=="continuous" & predictor.label!="ECCPD" & dep.var.type=="binary" ~ exp(b-1.96*SE)
      ,predictor.type=="continuous" & predictor.label=="ECCPD" & dep.var.type=="binary" ~ exp(caffeine.per.cup*(b-1.96*SE))
      ,predictor.type=="binary" & dep.var.type=="continuous" ~ log(2)*(b-1.96*SE)
      ,predictor.type=="binary" & dep.var.type=="binary" ~ exp(log(2)*(b-1.96*SE))
                                )
    ,effect.size.UBound=dplyr::case_when(
       predictor.type=="continuous" & predictor.label!="ECCPD" & dep.var.type=="continuous" ~ b+1.96*SE
      ,predictor.type=="continuous" & predictor.label=="ECCPD" & dep.var.type=="continuous" ~ caffeine.per.cup*(b+1.96*SE)
      ,predictor.type=="continuous" & predictor.label!="ECCPD" & dep.var.type=="binary" ~ exp(b+1.96*SE)
      ,predictor.type=="continuous" & predictor.label=="ECCPD" & dep.var.type=="binary" ~ exp(caffeine.per.cup*(b+1.96*SE))
      ,predictor.type=="binary" & dep.var.type=="continuous" ~ log(2)*(b+1.96*SE)
      ,predictor.type=="binary" & dep.var.type=="binary" ~ exp(log(2)*(b+1.96*SE))
                                 )
        ) # dim(coeff2) 373 16
    
# Customised sorting for the predictors
predictor.in.this.order <- c("Intercept"
                             ,"X20453_0_0_recoded","X20160_recode","caffeine.per.day","complete_alcohol_unitsweekly"
                             ,"age","inferred.sex","overall_health_rating","quali.edu.6138.recoded.z.score","TDI"
                             ,paste0("PC",c(1:10))) # "merged.pack.years.20161.z.scores" # length(predictor.in.this.order) 20

coeff2$predictor <- factor(coeff2$predictor
                           ,levels = predictor.in.this.order)

coeff2.sorted <- coeff2[order(coeff2$iteration,coeff2$dep.var.label,coeff2$predictor),] # dim(coeff2.sorted) 1362   16

# Give row number to dictate sorting
coeff2.sorted$row.order <- c(1:nrow(coeff2.sorted)) # dim(coeff2.sorted)

# Export merged phenotype data
ExportFileTabSeparated(data = ukb.phenos.covar.nonwhite.related.ID2.removed 
                       ,output.file.path = paste0(loc.obs.assoc,"phenotypes_combined.tsv")) # dim(ukb.phenos.covar.nonwhite.related.ID2.removed) 361283     32 

# Export regression results
ExportFileTabSeparated(data = coeff2.sorted
                       ,output.file.path = paste0(loc.obs.assoc,"phenotypic-association-results_outcomes-predictors-matched-MR-outcomes-exposures.tsv"))

# Export sample sizes
sample.sizes <- dplyr::bind_rows(base.logis.sample.size,base.linea.sample.size) # dim(sample.sizes) 21 6

ExportFileTabSeparated(data = sample.sizes
                       ,output.file.path = paste0(loc.obs.assoc,"sample-sizes_phenotypic-associations.tsv"))

#----------------------------------------------------------------------------------------------
#-------------------------------- This is the end of this program
#----------------------------------------------------------------------------------------------




#------------------------------------------------------------
# Why is p value shown as 0 for inferred.sex in the linear regression of ESDPW?
#------------------------------------------------------------
table(data$inferred.sex)
binary.var.names
contin.var.names

formula_i5 <- "complete_alcohol_unitsweekly~factor(X20453_0_0_recoded)+factor(X20160_recode)+caffeine.per.day+merged_pack_years_20161+age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8+ PC9+ PC10 + TDI + overall_health_rating + factor(inferred.sex) + quali.edu.6138.recoded.z.score"

test <- lm(data, formula = formula_i5)
summary(test)


#-------------------------- update the following code if Model 6 to 9 need to be redone-----------------------------#


# Set up parts of equations
equation.covariates <- "age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +TDI + overall_health_rating + factor(inferred.sex) + quali.edu.6138.recoded.z.score" # 

formula.1 <- paste0("X20160_recode ~ merged_pack_years_20161 + caffeine.per.day + complete_alcohol_unitsweekly + ",equation.covariates)
formula.2 <- paste0("X20160_recode ~ merged.pack.years.20161.z.scores + caffeine.per.day + complete_alcohol_unitsweekly + ",equation.covariates)

model.1 <- glm(formula = formula.1, data=data, family = binomial(link="logit"))
model.2 <- glm(formula = formula.2, data=data, family = binomial(link="logit"))

names(data)

table(data$X20160_recode)


pdf("test.pdf")
boxplot(merged_pack_years_20161~X20160_recode,data = data)
dev.off()
psych::describeBy(data$merged_pack_years_20161 , data$X20160_recode)


dd<-data[,c("merged_pack_years_20161","X20160_recode")]


summary(model.1)
# Get coefficients
model.1[["coefficients"]]
model.2[["coefficients"]] 




#----------------------------------------------------------
# Print results from selective columns for copy-and-paste
#----------------------------------------------------------
# Model 1: outcome=CI, predictors include other SU traits, covariates
formula.per.outcome[1,"formula"]
# Model 1 result
base.logis.coeff[,c(1,2,3,5)]

#-------------------------------------------------------------------------------------
# Model 5: relationship between cannabis initiation (as outcome) and caffeine consumed per day (as predictor) with 2 SNPs (rs2472297,rs4410790) added as covariates
#-------------------------------------------------------------------------------------
# Create a formula to construct the logistic model
outcome <- "X20453_0_0_recoded"

covariates <- "age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + TDI + overall_health_rating + factor(inferred.sex) + quali.edu.6138.recoded.z.score + factor(X20116_recodeFinal) + rs2472297 + rs4410790" # factor(EA.6138.uniCollegeDeg)

predictors <- c("caffeine.per.day+merged_pack_years_20161+complete_alcohol_unitsweekly")

formula <- paste0(outcome,"~", predictors,"+",covariates)

# Run logistic modelling
data <- ukb.phenos.covar.nonwhite.related.ID2.removed # dim(data) 361283     27

logi.fit.outc.CI.expo.CCPD.predictor.covar.2SNPs <- glm(formula= formula,data=data,family=binomial(link = "logit"))

# Get model summary
logi.summ.outc.CI.expo.CCPD.predictor.covar.2SNPs <- summary(logi.fit.outc.CI.expo.CCPD.predictor.covar.2SNPs) # str(logi.summ.outc.CI.expo.CCPD.predictor.covar.2SNPs) a list of 18

# Get coefficients
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs <- logi.summ.outc.CI.expo.CCPD.predictor.covar.2SNPs[["coefficients"]] # dim(mod.logi.summ.coeff) 18 4

# Calculate number of observations used in the modelling
numb.obs.raw.data <- nrow(data)
numb.obs.deleted <- length(logi.summ.outc.CI.expo.CCPD.predictor.covar.2SNPs[["na.action"]])
numb.obs.analysed <- numb.obs.raw.data-numb.obs.deleted

# Convert a matrix to a data.frame keeping its dimnames
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df <- plyr::adply(logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs,1,c)
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df$X1 <- as.character(logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df$X1) # dim(logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df) 20 5

# Rename columns
colnames(logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df) <- c("predictor","estimate","SE","z.value","p.value")

# Add name of the dependent variable
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df$dep.var <- outcome
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df$dep.var.label <- "CI"

# Add iteration number 
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df$iteration <- 5 # dim(logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df) 20 8

# Display results from selective columns
logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df[,c(1:5)]

#-------------------------------------------------------------------------------------
# Model 6: relationship between cannabis initiation (as outcome) and age + sex + 8 PCs + 2 SNPs (rs2472297,rs4410790) added as covariates
#-------------------------------------------------------------------------------------
# Create a formula to construct the logistic model
outcome <- "X20453_0_0_recoded"

covariates <- "age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + factor(inferred.sex) + rs2472297 + rs4410790"

formula <- paste0(outcome,"~", covariates)

# Run logistic modelling
data <- ukb.phenos.covar.nonwhite.related.ID2.removed # dim(data) 361283     27

logi.fit.outc.CI.covar.2SNPs <- glm(formula= formula,data=data,family=binomial(link = "logit"))

# Get model summary
logi.summ.outc.CI.covar.2SNPs <- summary(logi.fit.outc.CI.covar.2SNPs) # str(logi.summ.outc.CI.covar.2SNPs) a list of 18

# Get coefficients
logi.coef.outc.CI.covar.2SNPs <- logi.summ.outc.CI.covar.2SNPs[["coefficients"]] # dim(logi.coef.outc.CI.covar.2SNPs) 13 4

# Calculate number of observations used in the modelling
numb.obs.raw.data <- nrow(data)
numb.obs.deleted <- length(logi.summ.outc.CI.covar.2SNPs[["na.action"]])
numb.obs.analysed <- numb.obs.raw.data-numb.obs.deleted

# Convert a matrix to a data.frame keeping its dimnames
logi.coef.df.outc.CI.covar.2SNPs <- plyr::adply(logi.coef.outc.CI.covar.2SNPs,1,c)
logi.coef.df.outc.CI.covar.2SNPs$X1 <- as.character(logi.coef.df.outc.CI.covar.2SNPs$X1) # dim(logi.coef.df.outc.CI.covar.2SNPs) 13 5

# Rename columns
colnames(logi.coef.df.outc.CI.covar.2SNPs) <- c("predictor","estimate","SE","z.value","p.value")

# Add name of the dependent variable
logi.coef.df.outc.CI.covar.2SNPs$dep.var <- outcome
logi.coef.df.outc.CI.covar.2SNPs$dep.var.label <- "CI"

# Add iteration number 
logi.coef.df.outc.CI.covar.2SNPs$iteration <- 6 # dim(logi.coef.df.outc.CI.covar.2SNPs) 13 8

# Display results from selective columns
logi.coef.df.outc.CI.covar.2SNPs[,c(1:5)]

#-------------------------------------------------------------------------------------
# Model 7: relationship between caffeine consumed per day (as outcome) and cannabis initiation (as predictor) with 2 SNPs (rs2472297,rs4410790) added as covariates
#-------------------------------------------------------------------------------------
# Create a formula to construct the logistic model
outcome <- "caffeine.per.day"

covariates <- "age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + TDI + overall_health_rating + factor(inferred.sex) + quali.edu.6138.recoded.z.score + factor(X20116_recodeFinal)+ rs2472297 + rs4410790" # factor(EA.6138.uniCollegeDeg)

predictors <- "factor(X20453_0_0_recoded) + merged_pack_years_20161 + complete_alcohol_unitsweekly"

formula <- paste0(outcome,"~", predictors,"+",covariates)

# Run linear regression
data <- ukb.phenos.covar.nonwhite.related.ID2.removed # dim(data) 361283     27

linear.fit.outc.CCPD.predictors.covar.2SNPs <- lm(formula= formula,data=data)

# Get model summary
linear.summ.outc.CCPD.predictors.covar.2SNPs <- summary(linear.fit.outc.CCPD.predictors.covar.2SNPs) # str(linear.summ.outc.CCPD.predictors.covar.2SNPs) a list of 18

# Get coefficients
linear.coef.outc.CCPD.predictors.covar.2SNPs <- linear.summ.outc.CCPD.predictors.covar.2SNPs[["coefficients"]] # dim(linear.coef.outc.CCPD.predictors.covar.2SNPs) 20 4

# Convert a matrix to a data.frame keeping its dimnames
linear.coef.df.outc.CCPD.predictors.covar.2SNPs <- plyr::adply(linear.coef.outc.CCPD.predictors.covar.2SNPs,1,c)
linear.coef.df.outc.CCPD.predictors.covar.2SNPs$X1 <- as.character(linear.coef.df.outc.CCPD.predictors.covar.2SNPs$X1) # dim(linear.coef.df.outc.CCPD.predictors.covar.2SNPs) 20 5

# Rename columns
colnames(linear.coef.df.outc.CCPD.predictors.covar.2SNPs) <- c("predictor","estimate","SE","t.value","p.value")

# Add name of the dependent variable
linear.coef.df.outc.CCPD.predictors.covar.2SNPs$dep.var <- outcome
linear.coef.df.outc.CCPD.predictors.covar.2SNPs$dep.var.label <- "ECCPD"

# Add iteration number 
linear.coef.df.outc.CCPD.predictors.covar.2SNPs$iteration <- 7 # dim(linear.coef.df.outc.CCPD.predictors.covar.2SNPs) 20 8

#-------------------------------------------------------------------------------------
# Model 8: relationship between CCPD (outcome) and covariates
#-------------------------------------------------------------------------------------
# Create a formula to construct the logistic model
outcome <- "caffeine.per.day"

covariates <- "age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + factor(inferred.sex) + rs2472297 + rs4410790"

formula <- paste0(outcome,"~", covariates)

# Run linear regression
data <- ukb.phenos.covar.nonwhite.related.ID2.removed # dim(data) 361443     25

linear.fit.outc.CCPD.covar.2SNPs <- lm(formula= formula,data=data)

# Get model summary
linear.summ.outc.CCPD.covar.2SNPs <- summary(linear.fit.outc.CCPD.covar.2SNPs) # str(linear.summ.outc.CCPD.covar.2SNPs) a list of 18

# Get coefficients
linear.coef.outc.CCPD.covar.2SNPs <- linear.summ.outc.CCPD.covar.2SNPs[["coefficients"]] # dim(linear.coef.outc.CCPD.covar.2SNPs) 13 4

# Convert a matrix to a data.frame keeping its dimnames
linear.coef.df.outc.CCPD.covar.2SNPs <- plyr::adply(linear.coef.outc.CCPD.covar.2SNPs,1,c)
linear.coef.df.outc.CCPD.covar.2SNPs$X1 <- as.character(linear.coef.df.outc.CCPD.covar.2SNPs$X1) # dim(linear.coef.df.outc.CCPD.covar.2SNPs) 13 5

# Rename columns
colnames(linear.coef.df.outc.CCPD.covar.2SNPs) <- c("predictor","estimate","SE","t.value","p.value")
# Add name of the dependent variable
linear.coef.df.outc.CCPD.covar.2SNPs$dep.var <- outcome
linear.coef.df.outc.CCPD.covar.2SNPs$dep.var.label <- "ECCPD"

# Add iteration number 
linear.coef.df.outc.CCPD.covar.2SNPs$iteration <- 8 # dim(linear.coef.df.outc.CCPD.covar.2SNPs) 13 8

# Display results from selective columns
linear.coef.df.outc.CCPD.covar.2SNPs[,c(1:5)]

#-----------------------------------------------------------------------------------
# Combine all modelling results above 
#-----------------------------------------------------------------------------------

#----------------------------------------------
# Part 1: Combine model 5, 6, 7, 8
#----------------------------------------------

# Combine logistic modelling results (Model 5 and 6) 
# Calculate odds ratios and 95% CI for coefficients
postHoc.logi.mod.summ <- dplyr::bind_rows(logi.coef.outc.CI.expo.CCPD.predictor.covar.2SNPs.df
                                                   ,logi.coef.df.outc.CI.covar.2SNPs) %>% 
  dplyr::mutate(effect.size=exp(estimate)
                ,effect.size.lower.bound=exp(estimate - 1.96*SE)
                ,effect.size.upper.bound=exp(estimate + 1.96*SE)) # dim(postHoc.logi.mod.summ) 33 11

# Combine linear modelling results (Model 7 and 8) 
# Calculate odds ratios and 95% CI for coefficients
postHoc.linear.mod.summ <- dplyr::bind_rows(linear.coef.df.outc.CCPD.predictors.covar.2SNPs
                                                     ,linear.coef.df.outc.CCPD.covar.2SNPs) %>%
  dplyr::mutate( effect.size=estimate
                ,effect.size.lower.bound= estimate - 1.96*SE
                ,effect.size.upper.bound= estimate + 1.96*SE) # dim(postHoc.linear.mod.summ) 33 11

# Format values to 3 decimal places 
## Logistic regression results
tem1 <- postHoc.logi.mod.summ
tem1[,c("estimate","SE","effect.size","effect.size.lower.bound","effect.size.upper.bound")] <- lapply(tem1[,c("estimate","SE","effect.size","effect.size.lower.bound","effect.size.upper.bound")],Round.numeric.3decimal)

## Linear regression results
tem2 <- postHoc.linear.mod.summ
tem2[,c("estimate","SE","effect.size.lower.bound","effect.size.upper.bound")] <- lapply(tem2[,c("estimate","SE","effect.size.lower.bound","effect.size.upper.bound")],Round.numeric.3decimal)

# Add a commonly-named column with wanted estimates from both regression results
## For logistic, it is OR (95% CI)
## For linear, it is beta (95% CI)
tem1$effect.size.95CI <- with(tem1,paste0(effect.size," [",effect.size.lower.bound,", ",effect.size.upper.bound,"]"))
tem1$model <- "logistic regression"
tem1$effect.size.name <- "odds ratio (95% CI)" # dim(tem1) 33 14

tem2$effect.size.95CI <- with(tem2,paste0(estimate," [",effect.size.lower.bound,", ",effect.size.upper.bound,"]"))
tem2$model <- "linear regression"
tem2$effect.size.name <- "beta (95% CI)" # dim(tem2) 33 13

# Drop columns for merging the 2 results
tem1.small <- tem1 %>% dplyr::select(-c(z.value)) # dim(tem1.small) 33 13
tem2.small <- tem2 %>% dplyr::select(-c(t.value)) # dim(tem2.small) 33 13

# Combine the 2 regression results
## Delete 3 strings in a column (1) "factor(" , (2) ")1", (3) ")2"
strings.to.search <- glob2rx("^factor\\(|\\)1$|\\)2$")

tem1.tem2 <- dplyr::bind_rows(tem1.small,tem2.small) %>%
  dplyr::mutate(predictor= stringr::str_replace_all(predictor
                                                    ,pattern=strings.to.search
                                                    ,replacement = "")) # dim(tem1.tem2) 66 13

# Add labels for predictors by partial string matching
tem1.tem2 <- merge(x=tem1.tem2
                   ,y=variables_types[,c("varName","var.label")]
                   ,by.x="predictor"
                   ,by.y="varName"
                   ,all.x = TRUE)

tem1.tem2 <- data.table::setnames(tem1.tem2,old=c("var.label"),new=c("predictor.label")) # dim(tem1.tem2) 66 14

# Reorder columns
columns.order <- c("iteration","model","dep.var","dep.var.label","predictor","predictor.label","estimate","SE","effect.size.name","effect.size","effect.size.lower.bound","effect.size.upper.bound","effect.size.95CI","p.value") # length(columns.order) 14

tem1.tem2 <- tem1.tem2 %>% dplyr::select_(.dots = columns.order)

# Sort data 
tem1.tem2.sorted <- tem1.tem2[order(tem1.tem2$iteration,tem1.tem2$dep.var.label,tem1.tem2$predictor.label),] # dim(tem1.tem2.sorted) 66 14

#----------------------------------------------
# Part 2: rbind all results and export data
#----------------------------------------------

manu4.obs.assoc.results <- rbind(t1.t2.sorted,tem1.tem2.sorted) # dim(manu4.obs.assoc.results) 138 14




setwd(locScripts)
#file.copy("MR_step07_observational-analysis_effect-exposure-trait-on-outcome-trait_UKB-phenotype-data.R","zMR_step07_observational-analysis_effect-exposure-trait-on-outcome-trait_UKB-phenotype-data.R")

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#