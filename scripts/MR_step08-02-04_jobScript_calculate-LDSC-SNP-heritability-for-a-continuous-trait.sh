#!/bin/bash
## File path: ${locScripts}/MR_step08-02-04_jobScript_calculate-LDSC-SNP-heritability-for-a-continuous-trait.sh
## Modified from: ${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh
## Date created: 20190412
## Purpose: Calculate SNP heritability for a continuous trait. 
## Run dependency:
## Reference: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
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
file_path_munged=${v_file_path_munged}
consortium=${v_consortium}
trait=${v_trait}
folderPath_output=${v_folderPath_output}
# MF_1_population_prevalence=${v_MF_1_population_prevalence}
# MF_1_sample_prevalence=${v_MF_1_sample_prevalence}
# MF_2_filePath=${v_MF_2_filePath}
# MF_2_consortium=${v_MF_2_consortium}
# MF_2_trait=${v_MF_2_trait}
# MF_2_population_prevalence=${v_MF_2_population_prevalence}
# MF_2_sample_prevalence=${v_MF_2_sample_prevalence}

homeDir="/mnt/backedup/home/lunC";
locScripts=${homeDir}/scripts/MR_ICC_GSCAN_201806

#-------------------------------------------------------------------------------------------#
#---------Calculate SNP heritability in observed scale using LDSC---------------------------#
#-------------------------------------------------------------------------------------------#
module load ldsc

ldsc.py	--h2 ${file_path_munged} \
		--ref-ld-chr ${ref_ld_chr} \
		--w-ld-chr ${ref_ld_chr} \
		--out ${folderPath_output}/SNP-heritability_${consortium}-${trait}

#cp -n ${locScripts}/MR_step08-02-04_jobScript_calculate-LDSC-SNP-heritability-for-a-continuous-trait.sh ${locScripts}/MR_step08-02-05_jobScript_calculate-LDSC-SNP-heritability-for-a-binary-trait.sh
#---------------------------------------------------------------------------------------------------------------
##---------------------------------This is the end of this file-------------------------------------------------##
#---------------------------------------------------------------------------------------------------------------