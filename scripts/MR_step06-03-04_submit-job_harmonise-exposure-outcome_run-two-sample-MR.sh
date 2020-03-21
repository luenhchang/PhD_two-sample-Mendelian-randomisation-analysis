#!/bin/bash

# --------------------------------------------------------------------------------------------------------------------------------
# Program       : MR_step06-03-04_submit-job_harmonise-exposure-outcome_run-two-sample-MR.sh
# Modified from : /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step21-05-04_jobSubmit_2VarACE_genetic-corr-between-SUD-and-SUD-QIMR-adults.sh
# Date created  : 20190407
# Purpose       : Run a bash script containing commands to run R script (the job script) 
# How to run	: . ${locScripts}/MR_step06-03-04_submit-job_harmonise-exposure-outcome_run-two-sample-MR.sh > ${locHistory}/jobSubmitted_20190830_run_horizontal-pleiotropy-test
# Note			: The output of echo commands in the qsub script is sent to the psb out files
#----------------------------------------------------------------------------------------------------------------------------------
# Run dependency    : ${locScripts}/MR_step06-03-03_run-R-script-via-bash.sh
#----------------------------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input ${loc_input}/file-info_exposure-clumped-GWASs.tsv
# Input ${loc_input}/file-info_outcome-QCed-GWASs.tsv

# Outpu ${loc_output}/harmonised-data_exposure-clumped-GWAS-*-LDWindow-kb-10000-R2-0.01-p1-*-p2-*_outcome-* (180 files)
# Outpu ${loc_MR_result}/MR-analysis_exposure-*_outcome-*.tsv (180 files)
# Outpu ${loc_MR_report}/TwoSampleMR.Exposure*_against_Outcome*.html (210 files)
# Outpu ${loc_pleiotropy}/MR-Egger-intercept_exposure-*_outcome-*.tsv (210 files)
#--------------------------------------------------------------------------------------------------------------

# Sys.Date()History
#----------------------------------------------------------------------------------------------------------------------------------
# 20190830 qsub jobs 8772054-8772273 (210 jobs) Run only horizontal pleiotropy tests
# 20190823 qsub jobs 8710496-8710715 (210 jobs) Generated only MR report HTML files
# 20190812 qsub jobs 8604218-8604438 (221 jobs)
# 20190408 qsub jobs 7307910-7308089 (180 jobs). Got exposure and outcome harmonised. Done MR analysis on harmonised data
# 20190408 qsub jobs 7302854-7303033 (180 jobs)
#----------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------
# Folder locations under my home diretory
#-----------------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/MR_step06-03-03_run-R-script-via-bash.sh";

#-----------------------------------------------------
# Folder locations under my working diretory
#-----------------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806";
loc_input="${locMR}/two-sample-MR/input" # Location of input files
# Location of output files
loc_output="${loc_input}/harmonised-data"
loc_MR_result=${locMR}/two-sample-MR/output
loc_MR_result_archive=${loc_MR_result}/archive
pbs_output_dir="${loc_output}/pbs_output"
# This folder is where all the MR report HTML files are exported
loc_MR_report=${locMR}/two-sample-MR/output_MR-reports;

# This folder is where all the horizontal pleiotropy test results are exported
loc_pleiotropy=${locMR}/two-sample-MR/output_horizontal-pleiotropy;

mkdir -p ${loc_output} ${pbs_output_dir} ${loc_MR_result} ${loc_MR_result_archive} ${loc_MR_report} ${loc_pleiotropy};

echo ${homeDir}
echo ${locScripts}
echo ${locHistory}
echo ${jobScriptFilePath}

#------------------------------------------------------------------------
#---Move existing output files to an archive folder before a reanalysis
## Comment out the code after use
#------------------------------------------------------------------------
#find ${loc_MR_result} -name "MR-analysis_exposure-*.tsv" -exec mv {} ${loc_MR_result_archive} \;

#-------------------------------------------------------
# Set up resources requested for submitting PBS jobs
#-------------------------------------------------------
num_cpu=1;
#runTime_requested=2:00:00;
runTime_requested=1:00:00;  
memory_requested=5gb;
node_brand_name=Intel; # This avoids node 34

# Location of 2 TSV files to loop thru each line of them, excluding their header rows
filePath_tsv_1="${loc_input}/file-info_exposure-clumped-GWASs.tsv";
filePath_tsv_2="${loc_input}/file-info_outcome-QCed-GWASs.tsv";

#-------------------------------------------------------------------------------------------
# Loop through all combinations of $filePath_tsv_1 (22 lines) and $filePath_tsv_2 (11 lines)
#-------------------------------------------------------------------------------------------
## Number of total iterations: 242 (22*11)
### Number of iterations that are set to not run: 22 (11 traits* 2 p value thresholds for clumping)
### Number of iterations that are set to run: 221 (242-22)
IFS=$'\n';
count=0;
for lineF1 in `tail -n+2 $filePath_tsv_1`;do
	# Pass column values of file 1 to variables
	tsv_1_filePath=$(echo $lineF1 | cut -d$'\t' -f1);
	tsv_1_fileNameSpecial=$(echo $lineF1 | cut -d$'\t' -f2);
	tsv_1_consortium=$(echo $lineF1 | cut -d$'\t' -f3);
	tsv_1_trait=$(echo $lineF1 | cut -d$'\t' -f4);
	tsv_1_file_delimiter=$(echo $lineF1 | cut -d$'\t' -f5);
	tsv_1_colname_for_SNP=$(echo $lineF1 | cut -d$'\t' -f6);
	tsv_1_colname_for_beta=$(echo $lineF1 | cut -d$'\t' -f7);
	tsv_1_colname_for_SE=$(echo $lineF1 | cut -d$'\t' -f8);
	tsv_1_colname_for_effect_allele=$(echo $lineF1 | cut -d$'\t' -f9);
	tsv_1_colname_for_other_allele=$(echo $lineF1 | cut -d$'\t' -f10);
	tsv_1_colname_for_P_value=$(echo $lineF1 | cut -d$'\t' -f11);
	tsv_1_clumping_p1_value=$(echo $lineF1 | cut -d$'\t' -f12);
	tsv_1_clumping_p2_value=$(echo $lineF1 | cut -d$'\t' -f13);
	
	for lineF2 in `tail -n+2 $filePath_tsv_2`;do
		count=$((${count}+1));
		jobName="harmonise_exposure_outcome_pair${count}";
		# Pass column values of file 2 to variables
		tsv_2_filePath=$(echo $lineF2 | cut -d$'\t' -f1);
		tsv_2_consortium=$(echo $lineF2 | cut -d$'\t' -f2);
		tsv_2_trait=$(echo $lineF2 | cut -d$'\t' -f3);
		tsv_2_file_delimiter=$(echo $lineF2 | cut -d$'\t' -f4);
		tsv_2_colname_for_SNP=$(echo $lineF2 | cut -d$'\t' -f5);
		tsv_2_colname_for_beta=$(echo $lineF2 | cut -d$'\t' -f6);
		tsv_2_colname_for_SE=$(echo $lineF2 | cut -d$'\t' -f7);
		tsv_2_colname_for_effect_allele=$(echo $lineF2 | cut -d$'\t' -f8);
		tsv_2_colname_for_other_allele=$(echo $lineF2 | cut -d$'\t' -f9);
		tsv_2_colname_for_P_value=$(echo $lineF2 | cut -d$'\t' -f10);
		# Inspect variables passed from file 1
		echo "=============================================== iteration ${count} ===============================";
		echo "tsv_1_filePath=$tsv_1_filePath"; 
		echo "tsv_1_fileNameSpecial=$tsv_1_fileNameSpecial"
		echo "tsv_1_consortium=$tsv_1_consortium";
		echo "tsv_1_trait=$tsv_1_trait";
		echo "tsv_1_file_delimiter=$tsv_1_file_delimiter";
		echo "tsv_1_colname_for_SNP=$tsv_1_colname_for_SNP";
		echo "tsv_1_colname_for_beta=$tsv_1_colname_for_beta";
		echo "tsv_1_colname_for_SE=$tsv_1_colname_for_SE";
		echo "tsv_1_colname_for_effect_allele=$tsv_1_colname_for_effect_allele";
		echo "tsv_1_colname_for_other_allele=$tsv_1_colname_for_other_allele";
		echo "tsv_1_colname_for_P_value=$tsv_1_colname_for_P_value";
		echo "tsv_1_clumping_p1_value=$tsv_1_clumping_p1_value";
		echo "tsv_1_clumping_p2_value=$tsv_1_clumping_p2_value";
		# Inspect these variables passed from file 2
		echo "tsv_2_filePath=$tsv_2_filePath"; 
		echo "tsv_2_consortium=$tsv_2_consortium";
		echo "tsv_2_trait=$tsv_2_trait";
		echo "tsv_2_file_delimiter=$tsv_2_file_delimiter"
		echo "tsv_2_colname_for_SNP=$tsv_2_colname_for_SNP";
		echo "tsv_2_colname_for_beta=$tsv_2_colname_for_beta";
		echo "tsv_2_colname_for_SE=$tsv_2_colname_for_SE";
		echo "tsv_2_colname_for_effect_allele=$tsv_2_colname_for_effect_allele";
		echo "tsv_2_colname_for_other_allele=$tsv_2_colname_for_other_allele";
		echo "tsv_2_colname_for_P_value=$tsv_2_colname_for_P_value";
		# Create output file paths
		filePath_output_harmonised_data="${loc_output}/harmonised-data_exposure-${tsv_1_fileNameSpecial}_outcome-${tsv_2_consortium}-${tsv_2_trait}.tsv";
		filePath_output_MR_analysis="${loc_MR_result}/MR-analysis_exposure-${tsv_1_consortium}-${tsv_1_trait}-${tsv_1_clumping_p1_value}-${tsv_1_clumping_p2_value}_outcome-${tsv_2_consortium}-${tsv_2_trait}.tsv"
		filePath_output_hori_pleio="${loc_pleiotropy}/MR-Egger-intercept_exposure-${tsv_1_consortium}-${tsv_1_trait}-${tsv_1_clumping_p1_value}-${tsv_1_clumping_p2_value}_outcome-${tsv_2_consortium}-${tsv_2_trait}.tsv"

		echo "filePath_output_harmonised_data=${filePath_output_harmonised_data}";
		echo "filePath_output_MR_analysis=$filePath_output_MR_analysis" ;
		echo "filePath_output_hori_pleio=${filePath_output_hori_pleio}" ;

		# Conditionally qsub the shell script that runs a Rscript
		if [ "$tsv_1_consortium" = "$tsv_2_consortium" ] && [ "$tsv_1_trait" = "$tsv_2_trait" ] ;
			then echo "No action taken on identical exposure and outcome";
		else 
		# Test qsub command for one iteration
		echo "qsub -N $jobName -v v_tsv_1_filePath=${tsv_1_filePath},v_tsv_1_consortium=${tsv_1_consortium},v_tsv_1_trait=${tsv_1_trait},v_tsv_1_file_delimiter=${tsv_1_file_delimiter},v_tsv_1_colname_for_SNP=${tsv_1_colname_for_SNP},v_tsv_1_colname_for_beta=${tsv_1_colname_for_beta},v_tsv_1_colname_for_SE=${tsv_1_colname_for_SE},v_tsv_1_colname_for_effect_allele=${tsv_1_colname_for_effect_allele},v_tsv_1_colname_for_other_allele=${tsv_1_colname_for_other_allele},v_tsv_1_colname_for_P_value=${tsv_1_colname_for_P_value},v_tsv_1_clumping_p1_value=${tsv_1_clumping_p1_value},v_tsv_1_clumping_p2_value=${tsv_1_clumping_p2_value},v_tsv_2_filePath=${tsv_2_filePath},v_tsv_2_consortium=${tsv_2_consortium},v_tsv_2_trait=${tsv_2_trait},v_tsv_2_file_delimiter=${tsv_2_file_delimiter},v_tsv_2_colname_for_SNP=${tsv_2_colname_for_SNP},v_tsv_2_colname_for_beta=${tsv_2_colname_for_beta},v_tsv_2_colname_for_SE=${tsv_2_colname_for_SE},v_tsv_2_colname_for_effect_allele=${tsv_2_colname_for_effect_allele},v_tsv_2_colname_for_other_allele=${tsv_2_colname_for_other_allele},v_tsv_2_colname_for_P_value=${tsv_2_colname_for_P_value},v_iteration=${count},v_filePath_output_harmonised_data=${filePath_output_harmonised_data},v_filePath_output_MR_analysis=${filePath_output_MR_analysis},v_folderPath_output_MR_report=${loc_MR_report},v_filePath_output_hori_pleio=${filePath_output_hori_pleio} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${node_brand_name} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
		# Run a R script multiple times, each as a batch job
		qsub -N $jobName -v v_tsv_1_filePath=${tsv_1_filePath},v_tsv_1_consortium=${tsv_1_consortium},v_tsv_1_trait=${tsv_1_trait},v_tsv_1_file_delimiter=${tsv_1_file_delimiter},v_tsv_1_colname_for_SNP=${tsv_1_colname_for_SNP},v_tsv_1_colname_for_beta=${tsv_1_colname_for_beta},v_tsv_1_colname_for_SE=${tsv_1_colname_for_SE},v_tsv_1_colname_for_effect_allele=${tsv_1_colname_for_effect_allele},v_tsv_1_colname_for_other_allele=${tsv_1_colname_for_other_allele},v_tsv_1_colname_for_P_value=${tsv_1_colname_for_P_value},v_tsv_1_clumping_p1_value=${tsv_1_clumping_p1_value},v_tsv_1_clumping_p2_value=${tsv_1_clumping_p2_value},v_tsv_2_filePath=${tsv_2_filePath},v_tsv_2_consortium=${tsv_2_consortium},v_tsv_2_trait=${tsv_2_trait},v_tsv_2_file_delimiter=${tsv_2_file_delimiter},v_tsv_2_colname_for_SNP=${tsv_2_colname_for_SNP},v_tsv_2_colname_for_beta=${tsv_2_colname_for_beta},v_tsv_2_colname_for_SE=${tsv_2_colname_for_SE},v_tsv_2_colname_for_effect_allele=${tsv_2_colname_for_effect_allele},v_tsv_2_colname_for_other_allele=${tsv_2_colname_for_other_allele},v_tsv_2_colname_for_P_value=${tsv_2_colname_for_P_value},v_iteration=${count},v_filePath_output_harmonised_data=${filePath_output_harmonised_data},v_filePath_output_MR_analysis=${filePath_output_MR_analysis},v_folderPath_output_MR_report=${loc_MR_report},v_filePath_output_hori_pleio=${filePath_output_hori_pleio} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${node_brand_name} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};
		fi
	done
done

#cp -n /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-04_submit-job_harmonise-exposure-outcome_run-two-sample-MR.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-05-04_submit-job_run-MR-PRESSO.sh  
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#