#!/bin/bash

# --------------------------------------------------------------------------------------------------------------------------------
# Program       : MR_step06-05-04_submit-job_run-MR-PRESSO.sh
# Modified from : /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-04_submit-job_harmonise-exposure-outcome_run-two-sample-MR.sh
# Date created  : 20190419
# Purpose       : Run a bash script containing commands to run R script (the job script) 
# How to run	: . ${locScripts}/MR_step06-05-04_submit-job_run-MR-PRESSO.sh > ${locHistory}/jobSubmitted_20191011_run-MR-PRESSO
# Note			: The output of echo commands in the qsub script is sent to the psb out files
#----------------------------------------------------------------------------------------------------------------------------------
# Run dependency    : ${locScripts}/MR_step06-05-03_run-R-script-via-bash.sh
#----------------------------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input ${loc_input}/file-info_harmonised-data-selected-pairs-exposures-outcomes.tsv
# Outpu ${loc_output}/MRPRESSO-global-test_exposure*_outcome-ICC-CI.tsv (4 files)
# Outpu ${loc_output}/MRPRESSO-outlier-test_exposure*_outcome-ICC-CI.tsv (4 files)
# Outpu ${loc_output}/MRPRESSO-distortion-test_exposure*_outcome-ICC-CI.tsv (4 files)
#--------------------------------------------------------------------------------------------------------------

# Sys.Date()History
#----------------------------------------------------------------------------------------------------------------------------------
# 20191011 qsub jobs 9302368-9302371 (4 jobs)	
# 20190822 qsub jobs 8703412-8703416 (5 jobs)
# 20190818 qsub jobs 8670988-8670991 (4 jobs)
# 20190714 qsub jobs 8367560 (1 job UKB-CCPD 1e-5 on ICC-CI)
# 20190420 qsub jobs 7487906-7487908 (3 jobs)
# 20190419 qsub jobs 7487405-7487407 (3 jobs)
#----------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";
jobScriptFilePath="${locScripts}/MR_step06-05-03_run-R-script-via-bash.sh";

#-------------------------------------------
# Folder locations under my working directory
#-------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806";
# Location of input files
loc_input="${locMR}/MR-PRESSO/input"
# Location of output files
loc_output=${locMR}/MR-PRESSO/output/
loc_output_archive=${loc_output}/archive
pbs_output_dir="${loc_output}/pbs_output"
mkdir -p ${pbs_output_dir};

#---------------------------------------------------
# Archive existing output files
#---------------------------------------------------
# ls ${loc_output}/MRPRESSO*.tsv | wc -l # 12
#find ${loc_output} -name "MRPRESSO-*.tsv" -exec mv {} ${loc_output_archive} \;

#---------------------------------------------------
# Set up resources requested for submitting PBS jobs
#---------------------------------------------------
num_cpu=1;
runTime_requested=1:00:00; 
memory_requested=2gb;
chip_brand=Intel;

# Location of TSV files to loop thru each line of them, excluding their header rows
filePath_tsv_1="${loc_input}/file-info_harmonised-data-selected-pairs-exposures-outcomes.tsv";

#-----------------------------------------------------------------------------
# Copy variable assignment from previous step to assign variables at this step
#-----------------------------------------------------------------------------
# file_path_harmonised_data=${v_file_path_harmonised_data}
# exposure_trait=${v_exposure_trait}
# outcome_trait=${v_outcome_trait}
# colname_SNP=${v_colname_SNP}
# colname_beta_outcome=${v_colname_beta_outcome}
# colname_beta_exposure=${v_colname_beta_exposure}
# colname_SD_outcome=${v_colname_SD_outcome}
# colname_SD_exposure=${v_colname_SD_exposure}
# numb_simulation=${v_numb_simulation}
# significance_threshold=${v_significance_threshold}
# clumping_p1=${v_clumping_p1}
# output_folder_path=${v_output_folder_path}

#--------------------------------------------------------------------
# Loop through each line of $filePath_tsv_1
#--------------------------------------------------------------------
## Number of iterations: 5
IFS=$'\n';
count=0;
for line in `tail -n+2 $filePath_tsv_1`;do
	# Pass column values of tsv file 1 to variables
	file_path_harmonised_data=$(echo $line | cut -d$'\t' -f1);
	exposure_trait=$(echo $line | cut -d$'\t' -f2);
	clumping_p1=$(echo $line | cut -d$'\t' -f3);
	outcome_trait=$(echo $line | cut -d$'\t' -f4);
	colname_SNP=$(echo $line | cut -d$'\t' -f5);
	colname_beta_outcome=$(echo $line | cut -d$'\t' -f6);
	colname_beta_exposure=$(echo $line | cut -d$'\t' -f7);
	colname_SD_outcome=$(echo $line | cut -d$'\t' -f8);
	colname_SD_exposure=$(echo $line | cut -d$'\t' -f9);
	numb_simulation=$(echo $line | cut -d$'\t' -f10);
	significance_threshold=$(echo $line | cut -d$'\t' -f11);
	count=$((${count}+1));
	jobName="MRPRESSO_exposure_outcome_pair${count}";
	# Inspect variables passed from file 1
	echo "=============================================== iteration ${count} ===============================";
	echo "file_path_harmonised_data=$file_path_harmonised_data"; 
	echo "exposure_trait=$exposure_trait"
	echo "clumping_p1=$clumping_p1";
	echo "outcome_trait=$outcome_trait";
	echo "colname_SNP=$colname_SNP";
	echo "colname_beta_outcome=$colname_beta_outcome";
	echo "colname_beta_exposure=$colname_beta_exposure";
	echo "colname_SD_outcome=$colname_SD_outcome";
	echo "colname_SD_exposure=$colname_SD_exposure";
	echo "numb_simulation=$numb_simulation";
	echo "significance_threshold=$significance_threshold";
	# Specify folder path for output file in R
	output_folder_path=${loc_output}
	# Test qsub command for one iteration
	echo "qsub -N $jobName -v v_file_path_harmonised_data=${file_path_harmonised_data},v_exposure_trait=${exposure_trait},v_clumping_p1=${clumping_p1},v_outcome_trait=${outcome_trait},v_colname_SNP=${colname_SNP},v_colname_beta_outcome=${colname_beta_outcome},v_colname_beta_exposure=${colname_beta_exposure},v_colname_SD_outcome=${colname_SD_outcome},v_colname_SD_exposure=${colname_SD_exposure},v_numb_simulation=${numb_simulation},v_significance_threshold=${significance_threshold},v_output_folder_path=${output_folder_path} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${chip_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath}";
	# Run a R script multiple times, each as a batch job
	qsub -N $jobName -v v_file_path_harmonised_data=${file_path_harmonised_data},v_exposure_trait=${exposure_trait},v_clumping_p1=${clumping_p1},v_outcome_trait=${outcome_trait},v_colname_SNP=${colname_SNP},v_colname_beta_outcome=${colname_beta_outcome},v_colname_beta_exposure=${colname_beta_exposure},v_colname_SD_outcome=${colname_SD_outcome},v_colname_SD_exposure=${colname_SD_exposure},v_numb_simulation=${numb_simulation},v_significance_threshold=${significance_threshold},v_output_folder_path=${output_folder_path} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${chip_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${jobScriptFilePath};
done

#cp -n /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-03-04_submit-job_harmonise-exposure-outcome_run-two-sample-MR.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step06-05-04_submit-job_run-MR-PRESSO.sh 
#----------------------------------------------------------------------------#
#----------------This is the end of this file--------------------------------#
#----------------------------------------------------------------------------#