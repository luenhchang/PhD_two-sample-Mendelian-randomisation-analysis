#!/bin/bash

# ---------------------------------------------------------------------------------------
# Program       : MR_step06-03-03_run-R-script-via-bash.sh
# Modified from : /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step21-05-03_run-R-script-via-bash.sh
# Date created  : 20190407
# Purpose       : Run a R scirpt (MR_step06-03-02) as a PBS job (submitted by MR_step06-03-04)
# How to run	: 
#----------------------------------------------------------------------------------------
# Run dependency  : MR_step06-03-02_jobScript-harmonise-exposure-outcome.R					
#------------------------------------------------------------------------------------------------------
# Sys.Date() History
#----------------------------------------------------------------------------------------------------------------------------
# 20190829	Because module load R/3.5.1 is used here, simply use library(package) in the downstream R scripts. Don't assign lib.loc
## Wrong: 	library("dplyr",lib.loc="/software/R/R-3.4.1/lib64/R/library")
## Correct: library("dplyr")
# 20190829	Added ${filePath_output_hori_pleio} argument for horizontal pleiotropy test results. Tests done using TwoSampleMR::mr_pleiotropy_test(harmonised.data)
# 20190823	Added a folder path for exporting MR report HTML files , generated using TwoSampleMR::mr_report(harmonised.data)	 
#-----------------------------------------------------------------------------------------------------------------------------

# Pass shell variables that are specified by qsub at MR_step06-03-04_submit-job_harmonise-exposure-outcome.sh to shell variables here
## On the right of equal sign: qsub variables from MR_step06-03-04_submit-job_harmonise-exposure-outcome.sh
## On the left of equal sign: Shell variables to pass to R script through Rscript command

tsv_1_filePath=${v_tsv_1_filePath}
tsv_1_consortium=${v_tsv_1_consortium}
tsv_1_trait=${v_tsv_1_trait}
tsv_1_file_delimiter=${v_tsv_1_file_delimiter}
tsv_1_colname_for_SNP=${v_tsv_1_colname_for_SNP}
tsv_1_colname_for_beta=${v_tsv_1_colname_for_beta}
tsv_1_colname_for_SE=${v_tsv_1_colname_for_SE}
tsv_1_colname_for_effect_allele=${v_tsv_1_colname_for_effect_allele}
tsv_1_colname_for_other_allele=${v_tsv_1_colname_for_other_allele}
tsv_1_colname_for_P_value=${v_tsv_1_colname_for_P_value}
tsv_1_clumping_p1_value=${v_tsv_1_clumping_p1_value}
tsv_1_clumping_p2_value=${v_tsv_1_clumping_p2_value}

tsv_2_filePath=${v_tsv_2_filePath}
tsv_2_consortium=${v_tsv_2_consortium}
tsv_2_trait=${v_tsv_2_trait}
tsv_2_file_delimiter=${v_tsv_2_file_delimiter}
tsv_2_colname_for_SNP=${v_tsv_2_colname_for_SNP}
tsv_2_colname_for_beta=${v_tsv_2_colname_for_beta}
tsv_2_colname_for_SE=${v_tsv_2_colname_for_SE}
tsv_2_colname_for_effect_allele=${v_tsv_2_colname_for_effect_allele}
tsv_2_colname_for_other_allele=${v_tsv_2_colname_for_other_allele}
tsv_2_colname_for_P_value=${v_tsv_2_colname_for_P_value}

iteration=${v_iteration}
filePath_output_harmonised_data=${v_filePath_output_harmonised_data}
filePath_output_MR_analysis=${v_filePath_output_MR_analysis}
folderPath_output_MR_report=${v_folderPath_output_MR_report}
filePath_output_hori_pleio=${v_filePath_output_hori_pleio}

# Inspect the variables that are passed from qsub -v to here
echo "=============================================================================================================================="
echo "========================================Inspect Shell variables passed by qsub command========================================"
echo "=============================================================================================================================="
echo "tsv_1_filePath=${tsv_1_filePath}"
echo "tsv_1_consortium=$tsv_1_consortium"
echo "tsv_1_trait=$tsv_1_trait"
echo "tsv_1_file_delimiter=$tsv_1_file_delimiter"
echo "tsv_1_colname_for_SNP=$tsv_1_colname_for_SNP"
echo "tsv_1_colname_for_beta=$tsv_1_colname_for_beta"
echo "tsv_1_colname_for_SE=${tsv_1_colname_for_SE}"
echo "tsv_1_colname_for_effect_allele=${tsv_1_colname_for_effect_allele}"
echo "tsv_1_colname_for_other_allele=${tsv_1_colname_for_other_allele}"
echo "tsv_1_colname_for_P_value=${tsv_1_colname_for_P_value}"
echo "tsv_1_clumping_p1_value=$tsv_1_clumping_p1_value"
echo "tsv_1_clumping_p2_value=$tsv_1_clumping_p2_value"
echo "tsv_2_filePath=${tsv_2_filePath}"
echo "tsv_2_consortium=${tsv_2_consortium}"
echo "tsv_2_trait=${tsv_2_trait}"
echo "tsv_2_file_delimiter=$tsv_2_file_delimiter"
echo "tsv_2_colname_for_SNP=${tsv_2_colname_for_SNP}"
echo "tsv_2_colname_for_beta=${tsv_2_colname_for_beta}"
echo "tsv_2_colname_for_SE=${tsv_2_colname_for_SE}"
echo "tsv_2_colname_for_effect_allele=${tsv_2_colname_for_effect_allele}"
echo "tsv_2_colname_for_other_allele=${tsv_2_colname_for_other_allele}"
echo "tsv_2_colname_for_P_value=${tsv_2_colname_for_P_value}"
echo "iteration=${iteration}"
echo "filePath_output_harmonised_data=${filePath_output_harmonised_data}"
echo "filePath_output_MR_analysis=$filePath_output_MR_analysis"
echo "folderPath_output_MR_report=${folderPath_output_MR_report}"
echo "filePath_output_hori_pleio=${filePath_output_hori_pleio}"

# Set up directory
locScripts="/mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806"
RScriptFileName="MR_step06-03-02_jobScript-harmonise-exposure-outcome_run-two-sample-MR.R"
RScriptFilePath=${locScripts}/${RScriptFileName}

# Load software R in order to run a R file through the Rscript command
module load R/3.5.1

# Run a R script using Rscript command
## ${RScriptFilePath} : path of the R script file to run
## arguments that will be passed into the R script file: ${trait1_name} ~ ${iteration}
echo "==============================================================================================================================="
echo "================================Pass arguments to R ==========================================================================="
echo "================================Run a R script using Rscript command =========================================================="
echo "==============================================================================================================================="
echo "Rscript --vanilla ${RScriptFilePath} ${tsv_1_filePath} ${tsv_1_consortium} ${tsv_1_trait} ${tsv_1_file_delimiter} ${tsv_1_colname_for_SNP} ${tsv_1_colname_for_beta} ${tsv_1_colname_for_SE} ${tsv_1_colname_for_effect_allele} ${tsv_1_colname_for_other_allele} ${tsv_1_colname_for_P_value} ${tsv_1_clumping_p1_value} ${tsv_1_clumping_p2_value} ${tsv_2_filePath} ${tsv_2_consortium} ${tsv_2_trait} ${tsv_2_file_delimiter} ${tsv_2_colname_for_SNP} ${tsv_2_colname_for_beta} ${tsv_2_colname_for_SE} ${tsv_2_colname_for_effect_allele} ${tsv_2_colname_for_other_allele} ${tsv_2_colname_for_P_value} ${iteration} ${filePath_output_harmonised_data} ${filePath_output_MR_analysis} ${folderPath_output_MR_report} ${filePath_output_hori_pleio}"  

Rscript --vanilla ${RScriptFilePath} ${tsv_1_filePath} ${tsv_1_consortium} ${tsv_1_trait} ${tsv_1_file_delimiter} ${tsv_1_colname_for_SNP} ${tsv_1_colname_for_beta} ${tsv_1_colname_for_SE} ${tsv_1_colname_for_effect_allele} ${tsv_1_colname_for_other_allele} ${tsv_1_colname_for_P_value} ${tsv_1_clumping_p1_value} ${tsv_1_clumping_p2_value} ${tsv_2_filePath} ${tsv_2_consortium} ${tsv_2_trait} ${tsv_2_file_delimiter} ${tsv_2_colname_for_SNP} ${tsv_2_colname_for_beta} ${tsv_2_colname_for_SE} ${tsv_2_colname_for_effect_allele} ${tsv_2_colname_for_other_allele} ${tsv_2_colname_for_P_value} ${iteration} ${filePath_output_harmonised_data} ${filePath_output_MR_analysis} ${folderPath_output_MR_report} ${filePath_output_hori_pleio}

# Copy this file for similar job
#cp -n /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-03_run-R-script-via-bash.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-05-03_run-R-script-via-bash.sh
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#