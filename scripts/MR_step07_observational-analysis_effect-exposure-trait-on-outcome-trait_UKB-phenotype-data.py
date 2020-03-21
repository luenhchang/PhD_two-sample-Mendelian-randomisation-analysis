#!/usr/bin/env python

##################################################################################
# Filename: MR_step07_observational-analysis_effect-exposure-trait-on-outcome-trait_UKB-phenotype-data.py
# Programmer: Chang
# Purpose: Replicate MR_step07_observational-analysis_effect-exposure-trait-on-outcome-trait_UKB-phenotype-data.R
# Date created: 20200131
# Dependent: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step00-01_recode-UKB-phenotypes-for-running-GWAS.R
# How to run this script in iPython
# %run /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step07_observational-analysis_effect-exposure-trait-on-outcome-trait_UKB-phenotype-data.py
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
#-----------------------------------------------------------------------------------------
import pandas as pd
from dfply import *

#-----------------------------------------------------------------------------------------
# Folder locations under Stuart MacGregor's lab
#-----------------------------------------------------------------------------------------
stuartmaDir = "/reference/data/UKBB_500k/versions/lab_stuartma/"
stuartma_pheno = stuartmaDir + "pheno/"
stuartma_pheno_coffee = stuartma_pheno + "coffee/"
stuartma_pheno_caffeine = stuartma_pheno_coffee + "ukb_1488-teaIntake_1498-coffeeIntake_1508-coffeeType/"

#-----------------------------------------------------------------------------------------
# Folder locations under my home directory
#-----------------------------------------------------------------------------------------
homeDir = "/mnt/backedup/home/lunC/"
locRFunction = homeDir + "scripts/RFunctions/"
locScripts = homeDir + "scripts/MR_ICC_GSCAN_201806/"
dir_ukbPheno = homeDir + "data/UKBiobank_phenotype/"

#-----------------------------------------------------------------------------------------
# Folder locations under my working directory
#-----------------------------------------------------------------------------------------
workingDir = "/mnt/lustre/working/lab_nickm/lunC/"
locMR = workingDir + "MR_ICC_GSCAN_201806/" # location of outcome data
loc_obs_assoc = locMR + "observational-associations/"

#-------------------------------------------------------------------------------------#
# Import a list of exposures and outcomes used in MR analyses
# Match predictors and outcomes from the UKB to the MR exposures and outcomes
#-------------------------------------------------------------------------------------#
#MR_exposure_outcome_traits = pd.read_csv('/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/observational-associations/exposure-outcome-pairs-used-in-MR-analyses.tsv', sep='\t', engine='python')

#print(MR_exposure_outcome_traits.shape) # (23, 4)

#-------------------------------------------------------------------------------------#
#------File paths of phenotypes, covariates, non-white, relatives---------------------#
#--- Copy recoded values and definition under each file path
#-------------------------------------------------------------------------------------#
# Education attainment
file_path_ukb_pheno_edu_att = stuartma_pheno + "coffee/additional_coffee_lifestyle.pheno"

# Estimated standard drinks per week
file_path_ukb_pheno_ESDPW = stuartmaDir + "collab/alcohol_BCAC/total_alcohol_unitsweekly.combined.pheno"

# Caffeine consumed per day (mg)
file_path_ukb_pheno_CCPD = stuartma_pheno_caffeine + "daily.caffeine.consumption.thru.coffee.tea"

# Covariates
file_path_ukb_pheno_covariates = stuartma_pheno + "may18_baseline_extra_covariate"

# Ever smoked
file_path_ukb_pheno_20160 = dir_ukbPheno + "ukb20160_everSmoked/ukb20160.phenoUtility.recoded"

# Pack years of smoking
file_path_ukb_pheno_20161 = stuartma_pheno + "smoking/pack_years.pheno"

# Ever taken cannabis
file_path_ukb_pheno_20453 = dir_ukbPheno + "/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded"

# Age at starting smoking in current smokers
file_path_ukb_pheno_3436 = dir_ukbPheno + "/ukb3436_ageStartedSmokingInCurrentSmokers/ukb3436.phenoUtility.recoded"

# Cigarettes per day
file_path_ukb_pheno_3456 = dir_ukbPheno + "/ukb3456_numCigareDaily/ukb3456.phenoUtility"

# IDs of non-caucasian
file_path_ukb_non_white = stuartmaDir + "exclude_related/notWhite.id"

# Relatives of UKB participants
file_path_ukb_relatives = stuartmaDir + "relatives_working/alltrait_rela"

#-------------------------------------------------------------------------------------#
#-----------------Import UKB phenotype files -----------------------------------------#
#-------------------------------------------------------------------------------------#
ukb_pheno_edu_att = pd.read_csv(file_path_ukb_pheno_edu_att, sep='\s+', engine='python', usecols=["FID","IID","qualification_edu_6138"])
print(ukb_pheno_edu_att.shape) # (487409, 3)






