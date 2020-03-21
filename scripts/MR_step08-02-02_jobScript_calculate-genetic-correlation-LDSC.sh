#!/bin/bash
## File path: ${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh 
## Modified from: ${locScripts}/MR_step04-01-01_jobScript_LD-based-SNP-clumping.sh
## date created: 20190411
## purpose: calculate genetic correlation using LDSC python script
## Run dependency: 
## How to run this script: 
## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20190813	Replaced ${ldsc_dir}/ldsc.py with module load ldsc ldsc.py
##--------------------------------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under Stuart's lab
#-------------------------------------------
ref_snplist=/reference/data/UKBB_500k/versions/lab_stuartma/LDSC_res/1000HGP/w_hm3.snplist
ldsc_dir=/reference/data/UKBB_500k/versions/lab_stuartma/scripts/ldsc
ref_ld_chr=/reference/data/UKBB_500k/versions/lab_stuartma/LDSC_res/1000HGP/eur_w_ld_chr/

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts=${homeDir}/scripts/MR_ICC_GSCAN_201806
locHistory=${homeDir}/history

# Pass qsub variables to shell variables
MF_1_filePath=${v_MF_1_filePath}
MF_1_consortium=${v_MF_1_consortium}
MF_1_trait=${v_MF_1_trait}
MF_1_population_prevalence=${v_MF_1_population_prevalence}
MF_1_sample_prevalence=${v_MF_1_sample_prevalence}
MF_2_filePath=${v_MF_2_filePath}
MF_2_consortium=${v_MF_2_consortium}
MF_2_trait=${v_MF_2_trait}
MF_2_population_prevalence=${v_MF_2_population_prevalence}
MF_2_sample_prevalence=${v_MF_2_sample_prevalence}
folderPath_output=${v_folderPath_output}

homeDir="/mnt/backedup/home/lunC";
locScripts=${homeDir}/scripts/MR_ICC_GSCAN_201806

#-------------------------------------------------------------------------------------------#
#------------------------------calculate rG using module ldsc-------------------------------#
#-------------------------------------------------------------------------------------------#
module load ldsc

ldsc.py	--rg ${MF_1_filePath},${MF_2_filePath} \
		--ref-ld-chr ${ref_ld_chr} \
		--w-ld-chr ${ref_ld_chr} \
		--out ${folderPath_output}/genetic-correlation_between_${MF_1_consortium}-${MF_1_trait}_and_${MF_2_consortium}-${MF_2_trait} \
		--samp-prev ${MF_1_sample_prevalence},${MF_2_sample_prevalence} \
		--pop-prev ${MF_1_population_prevalence},${MF_2_population_prevalence}

#cp -n ${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh ${locScripts}/MR_step08-02-04_jobScript_calculate-LDSC-SNP-heritability-for-a-continuous-trait.sh  		

#---------------------------------------------------------------------------------------------------------------
##---------------------------------This is the end of this file-------------------------------------------------##
#---------------------------------------------------------------------------------------------------------------