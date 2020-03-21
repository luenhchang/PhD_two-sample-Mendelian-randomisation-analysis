#!/bin/bash

# ---------------------------------------------------------------------------------------
# Program       : MR_step06-05-03_run-R-script-via-bash.sh
# Modified from : /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-03_run-R-script-via-bash.sh
# Date created  : 20190419
# Purpose       : Run a R scirpt (MR_step06-05-02) as a PBS job (submitted by MR_step06-05-04)
# How to run	: 
#----------------------------------------------------------------------------------------
# Run dependency: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-05-02_job-script_run-MR-PRESSO.R
#------------------------------------------------------------------------------------------------------
# Sys.Date() History
#------------------------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Pass shell variables that are specified by qsub at MR_step06-05-04_submit-job_run-MR-PRESSO.sh to shell variables here
## On the left of = sign: Shell variables to pass to R script through Rscript command. 
### Copy these variables from R replacing the dots with underscores
## On the right of = sign: qsub variables from MR_step06-05-04_submit-job_harmonise-exposure-outcome.sh
### Simply prefix the left with "v_" 
#-----------------------------------------------------------------------------------------------------------------------
file_path_harmonised_data=${v_file_path_harmonised_data}
exposure_trait=${v_exposure_trait}
outcome_trait=${v_outcome_trait}
colname_SNP=${v_colname_SNP}
colname_beta_outcome=${v_colname_beta_outcome}
colname_beta_exposure=${v_colname_beta_exposure}
colname_SD_outcome=${v_colname_SD_outcome}
colname_SD_exposure=${v_colname_SD_exposure}
numb_simulation=${v_numb_simulation}
significance_threshold=${v_significance_threshold}
clumping_p1=${v_clumping_p1}
output_folder_path=${v_output_folder_path}
 
# Inspect the variables that are passed from qsub -v to here
echo "=============================================================================================================================="
echo "========================================Inspect Shell variables passed by qsub command========================================"
echo "=============================================================================================================================="
echo "file_path_harmonised_data=${file_path_harmonised_data}"
echo "exposure_trait=${exposure_trait}"
echo "outcome_trait=${outcome_trait}"
echo "colname_SNP=${colname_SNP}"
echo "colname_beta_outcome=${colname_beta_outcome}"
echo "colname_beta_exposure=${colname_beta_exposure}"
echo "colname_SD_outcome=${colname_SD_outcome}"
echo "colname_SD_exposure=${colname_SD_exposure}"
echo "numb_simulation=${numb_simulation}"
echo "significance_threshold=${significance_threshold}"
echo "clumping_p1=${clumping_p1}"
echo "output_folder_path=${output_folder_path}"

#---------------------------------------------------------------------
# Set up directory
#---------------------------------------------------------------------
locScripts="/mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806"
#RScriptFileName="MR_step06-03-02_jobScript-harmonise-exposure-outcome_run-two-sample-MR.R"
RScriptFileName="MR_step06-05-02_job-script_run-MR-PRESSO.R"
RScriptFilePath=${locScripts}/${RScriptFileName}

# Load software R in order to run a R file through the Rscript command
module load R/3.4.1

# Run a R script using Rscript command
## ${RScriptFilePath} : path of the R script file to run
## arguments that will be passed into the R script file 
echo "==============================================================================================================================="
echo "================================Passing arguments to R ==========================================================================="
echo "================================Run a R script using Rscript command =========================================================="
echo "==============================================================================================================================="
echo "Rscript --vanilla ${RScriptFilePath} ${tsv_1_filePath} ${tsv_1_consortium} ${tsv_1_trait} ${tsv_1_file_delimiter} ${tsv_1_colname_for_SNP} ${tsv_1_colname_for_beta} ${tsv_1_colname_for_SE} ${tsv_1_colname_for_effect_allele} ${tsv_1_colname_for_other_allele} ${tsv_1_colname_for_P_value} ${tsv_1_clumping_p1_value} ${tsv_1_clumping_p2_value} ${tsv_2_filePath} ${tsv_2_consortium} ${tsv_2_trait} ${tsv_2_file_delimiter} ${tsv_2_colname_for_SNP} ${tsv_2_colname_for_beta} ${tsv_2_colname_for_SE} ${tsv_2_colname_for_effect_allele} ${tsv_2_colname_for_other_allele} ${tsv_2_colname_for_P_value} ${iteration} ${filePath_output_harmonised_data} ${filePath_output_MR_analysis}"  

Rscript --vanilla ${RScriptFilePath} ${file_path_harmonised_data} ${exposure_trait} ${outcome_trait} ${colname_SNP} ${colname_beta_outcome} ${colname_beta_exposure} ${colname_SD_outcome} ${colname_SD_exposure} ${numb_simulation} ${significance_threshold} ${clumping_p1} ${output_folder_path}

# Copy this file for similar job
#cp -n /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-03_run-R-script-via-bash.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-05-03_run-R-script-via-bash.sh
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#