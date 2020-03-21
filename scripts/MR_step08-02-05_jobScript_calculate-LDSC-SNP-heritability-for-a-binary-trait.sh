#!/bin/bash
## File path: ${locScripts}/MR_step08-02-05_jobScript_calculate-LDSC-SNP-heritability-for-a-binary-trait.sh
## Modified from: ${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh
## Date created: 20190412
## Purpose: Calculate SNP heritability for a binary trait 
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
population_prevalence=${v_population_prevalence}
sample_prevalence=${v_sample_prevalence}
folderPath_output=${v_folderPath_output}

#-------------------------------------------------------------------------------------------#
#---------Calculate SNP heritability in observed scale using LDSC---------------------------#
#-------------------------------------------------------------------------------------------#
module load ldsc

ldsc.py	--h2 ${file_path_munged} \
		--ref-ld-chr ${ref_ld_chr} \
		--w-ld-chr ${ref_ld_chr} \
		--samp-prev ${sample_prevalence} \
		--pop-prev ${population_prevalence} \
		--out ${folderPath_output}/SNP-heritability_${consortium}-${trait}

#cp -n ${locScripts}/MR_step08-02-05_jobScript_calculate-LDSC-SNP-heritability-for-a-binary-trait.sh
#---------------------------------------------------------------------------------------------------------------
##---------------------------------This is the end of this file-------------------------------------------------##
#---------------------------------------------------------------------------------------------------------------