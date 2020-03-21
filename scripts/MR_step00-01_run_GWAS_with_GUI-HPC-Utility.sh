#!/bin/bash

## File name: MR_step00-01_run_GWAS_with_GUI-HPC-Utility.sh
## Date created: 20180731
## Run dependency: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step00-01_recode-UKB-phenotypes-for-running-GWAS.R
## Note: (1) prior to this step, phenotype data need to be cleaned, recoded in MR_step00-01_recode-UKB-phenotypes-for-running-GWAS.R
## Purpose: (1) this script documents how GWAS is run using HPC_Utility.jar 

## Time 	Change
##--------------------------------------------------------------------------------------------------------
## 20191024 qsub jobs 9538938-9538961 (24 jobs for runnning GWAS on ever taken cannabis. Phenotype recoded as 2 and 1)
## 20191023	qsub jobs 9535541-9535564 (24 jobs for runnning GWAS on ever taken cannabis, No output generated because phenotype was coded as 1 and 0)
## 20190811	qsub jobs 8589162 (1 job for chr8. Chr8 was run on node34, a slow node. Avoid this node in future analysis)
## How to find out which node?
# qstat -fx 8587001|grep exec_host
#    exec_host = hpcnode034/4*8 
## 20190810	qsub jobs 8586994-8587015 (22 jobs, trait=caffeine.per.day)	
## 20180822	qsub jobs 6093306-6093327 (22 jobs, trait=all_coffee_cpd)
## 20180822	qsub jobs 6093284-6093305 (22 jobs, trait=merged_pack_years_20161) 	
## 20180820 qsub jobs 6087364-6087385 (22 jobs) 
## 20180731 qsub jobs 6038661-6038682 (22 jobs) 
## 20180731 qsub jobs 6038309-6038330 (22 jobs; an error coz phenotype file didn't begin with FID and then IID) 
##--------------------------------------------------------------------------------------------------------

# Output folders under /reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas/BOLT_LMM
outputMain="/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas";
outputFolderPath_ukb3456_NA20453="${outputMain}/BOLT_LMM/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis";
outputFolderPath_ESDPW_NA20453="${outputMain}/BOLT_LMM/UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis"; # folder created in PRS_UKB_201711_step00-00_recode_phenotype_for_GWAS.R
outputFolderPath_ukb20161_NA20453=${outputMain}/BOLT_LMM/UKB20161-pack-years-of-smoking_IID-NA-in-UKB20453-everUsedCannabis;
outputFolderPath_ukb_CCPD_NA20453=${outputMain}/BOLT_LMM/UKB-cups-of-coffee-per-day_IID-NA-in-UKB20453-everUsedCannabis;
output_folder_path_UKB_caffeine_NA20453=${outputMain}/BOLT_LMM/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis;
outputFolderPath_ukb20453=${outputMain}/plink2/UKB20453-ever-taken-cannabis;

phenotypeFolderPath_ukb3456_NA20453="${outputFolderPath_ukb3456_NA20453}/phenotype";

phenotypeFilePath_ESDPW_NA20453="${outputFolderPath_ESDPW_NA20453}/phenotype";
phenotypeFilePath_ukb20161_NA20453=${outputFolderPath_ukb20161_NA20453}/phenotype;
phenotypeFilePath_ukb_CCPD_NA20453=$outputFolderPath_ukb_CCPD_NA20453/phenotype;
phenotypeFilePath_ukb_caffeine_NA20453=${output_folder_path_UKB_caffeine_NA20453}/phenotype;

historyFolderPath_ukb3456_NA20453="${outputFolderPath_ukb3456_NA20453}/history";
historyFolderPath_ESDPW_NA20453="${outputFolderPath_ESDPW_NA20453}/history";
historyFilePath_ukb20161_NA20453=$outputFolderPath_ukb20161_NA20453/jobSubmitted_20180822_run-GWAS-UKB;
historyFilePath_ukb_CCPD_NA20453=$outputFolderPath_ukb_CCPD_NA20453/jobSubmitted_20180822_run-GWAS-UKB;
historyFilePath_ukb_caffeine_NA20453=${output_folder_path_UKB_caffeine_NA20453}/jobSubmitted_20190810_run-GWAS;
historyFilePath_ukb20453_everTakenCannabis=${outputFolderPath_ukb20453}/jobSubmitted_20191024_run-GWAS-binary-ever-taken-cannabis

mkdir -p ${outputFolderPath_ukb3456_NA20453} ${phenotypeFolderPath_ukb3456_NA20453} ${historyFolderPath_ukb3456_NA20453} ${outputFolderPath_ukb20453};

#----------------------------------------------------------------------------------------------------------------------
# Run HPC_Utility.jar to generate shell scripts for running GWAS
#----------------------------------------------------------------------------------------------------------------------
## step 1: in Windows cmd.exe go to the software folder typing "D:"
## step 2: type "cd D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_HPCUtility"
## step 3: run HPC_Utility.jar by typing "java -jar HPC_Utility.jar"
## step 4: a login window will pop up. Enter HPC username and password
## step 5: specify BOLT-LMM or plink2. Select output directory, phenotype directory, phenotype column name and click generate script
## step 6: run the script using the following command

## Document entry and output files from HPC Utility that will run GWAS analysis
## Type 	Def 			FilePath
##---------------------------------------------------------------------------------------------------------------------------
## Input	phenoFile		${phenotypeFolderPath_ukb3456_NA20453}/ukb3456_IID_NA_in_20453
## Input	phenoColname	X3456_mean 	 
## GUI		software		D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR_HPCUtility\HPC_Utility.jar
## Outpu 	script			${outputFolderPath_ukb3456_NA20453}/ukb3456_IID_NA_in_20453-BOLT-LMM_X3456_mean.sh
## Outpu 	script			${outputFolderPath_ukb3456_NA20453}/ukb3456_IID_NA_in_20453-BOLT-LMM_X3456_mean.prelim

## Input	phenoFile		${outputFolderPath_ESDPW_NA20453}/phenotype
## Input	phenoColname	complete_alcohol_unitsweekly
## Outpu 	script			${outputFolderPath_ESDPW_NA20453}/phenotype-BOLT-LMM_complete_alcohol_unitsweekly.sh
## Outpu 	script			${outputFolderPath_ESDPW_NA20453}/phenotype-BOLT-LMM_complete_alcohol_unitsweekly.prelim

## Input	phenoFile		$outputFolderPath_ukb20161_NA20453/phenotype
## Input	phenoColname	merged_pack_years_20161
## Outpu	script			${outputFolderPath_ukb20161_NA20453}/phenotype-BOLT-LMM_merged_pack_years_20161.sh
## Outpu	script			${outputFolderPath_ukb20161_NA20453}/phenotype-BOLT-LMM_merged_pack_years_20161.prelim

## Input	phenoFile		$outputFolderPath_ukb_CCPD_NA20453/phenotype
## Input	phenoColname	all_coffee_cpd
## Outpu	script			$outputFolderPath_ukb_CCPD_NA20453/phenotype-BOLT-LMM_all_coffee_cpd.sh
## Outpu	script			$outputFolderPath_ukb_CCPD_NA20453/phenotype-BOLT-LMM_all_coffee_cpd.prelim

## Input	phenoFile		${output_folder_path_UKB_caffeine_NA20453}/phenotype
## Input	phenoColname	caffeine.per.day
## Outpu	script			${output_folder_path_UKB_caffeine_NA20453}/phenotype-BOLT-LMM_caffeine.per.day.sh
## Outpu	script			${output_folder_path_UKB_caffeine_NA20453}/phenotype-BOLT-LMM_caffeine.per.day.prelim

## Input	phenoFile		/mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded
## Input	phenoColname	X20453_0_0_recoded
## Outpu	script			${outputFolderPath_ukb20453}/ukb20453.phenoUtility.recoded-plink2_X20453_0_0_recoded.sh
## Outpu	script			${outputFolderPath_ukb20453}/ukb20453.phenoUtility.recoded-plink2_X20453_0_0_recoded.prelim

## Copy GWAS on ever smoked to the plink2 folder

##--------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------
# Run GWAS on autosomes using scripts generated from HPC_Utility.jar 
#----------------------------------------------------------------------------------------------------------------------

. ${outputFolderPath_ukb3456_NA20453}/ukb3456_IID_NA_in_20453-BOLT-LMM_X3456_mean.sh > ${historyFolderPath_ukb3456_NA20453}/jobSubmitted_20180731_run_GWAS_BOLT-LMM_UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis

. ${outputFolderPath_ESDPW_NA20453}/phenotype-BOLT-LMM_complete_alcohol_unitsweekly.sh > ${historyFolderPath_ESDPW_NA20453}

. ${outputFolderPath_ukb20161_NA20453}/phenotype-BOLT-LMM_merged_pack_years_20161.sh > ${historyFilePath_ukb20161_NA20453}

. ${outputFolderPath_ukb_CCPD_NA20453}/phenotype-BOLT-LMM_all_coffee_cpd.sh > ${historyFilePath_ukb_CCPD_NA20453};

. ${output_folder_path_UKB_caffeine_NA20453}/phenotype-BOLT-LMM_caffeine.per.day.sh > ${historyFilePath_ukb_caffeine_NA20453}

. ${outputFolderPath_ukb20453}/ukb20453.phenoUtility.recoded-plink2_X20453_0_0_recoded.sh > ${historyFilePath_ukb20453_everTakenCannabis}

# Document output files after GWAS run by shell scripts submitted at previous step 
# GWAS file paths
#---------------------------------------------------------------------------------------------------------------------------
# ${outputFolderPath_ukb3456_NA20453}/BOLT-LMM-ukb3456_IID_NA_in_20453-pheno-X3456_mean/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_X3456_mean.bgen.assoc
# ${outputFolderPath_ESDPW_NA20453}/BOLT-LMM-phenotype-pheno-complete_alcohol_unitsweekly/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_complete_alcohol_unitsweekly.bgen.assoc

# ${outputFolderPath_ukb20161_NA20453}/BOLT-LMM-phenotype-pheno-merged_pack_years_20161/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_merged_pack_years_20161.bgen.assoc

# $outputFolderPath_ukb_CCPD_NA20453/BOLT-LMM-phenotype-pheno-all_coffee_cpd/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_all_coffee_cpd.bgen.assoc

# ${output_folder_path_UKB_caffeine_NA20453}/
#--------------------------------------------------------------------------------------------------------------------------------

