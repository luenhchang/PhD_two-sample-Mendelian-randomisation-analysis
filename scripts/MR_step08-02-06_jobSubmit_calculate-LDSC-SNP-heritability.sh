#!/bin/bash
 
## File name: ${locScripts}/MR_step08-02-06_jobSubmit_calculate-LDSC-SNP-heritability.sh
## Modified from: ${locScripts}/MR_step08-02-03_jobSubmit_calculate-genetic-correlation-LDSC.sh
## Date created: 20190412
## Note: 
## Purpose: calculate SNP heritability using LD score correlation python scripts
## Run dependency: 
### (1) ${locScripts}/MR_step08-02-04_jobScript_calculate-LDSC-SNP-heritability-for-a-continuous-trait.sh
### (2) ${locScripts}/MR_step08-02-05_jobScript_calculate-LDSC-SNP-heritability-for-a-binary-trait
### (3) module load ldsc ldsc.py

## How to run this file: . ${locScripts}/MR_step08-02-06_jobSubmit_calculate-LDSC-SNP-heritability.sh >  ${locHistory}/jobs-submitted_20190813_load-module-ldsc_calculate-SNP-heritability
## Note: Don't open the CSV files by right clicking edit with Excel. This may cause all columns appear in single column. Instead, import the CSV files in Excel locally, edit them and save them again as CSV, and transfer them back to remote folder. See instruction:  https://support.goteamup.com/hc/en-us/articles/203460681-FAQ-When-opening-CSV-file-all-data-appears-in-one-column

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input 	/reference/data/UKBB_500k/versions/lab_stuartma/LDSC_res/1000HGP/w_hm3.snplist	
# Input		${loc_LDSC_input}/file-info_munged-QCed-GWASs.tsv

# Outpu		${loc_LDSC_SNP_heritability}/SNP-heritability_*-*.log (11)
#------------------------------------------------------------------------------------------------------
## Time 	Change
##---------------------------------------------------------------------------------------------------------------------
## 20190813	qsub jobs 8616187-8616198 (11 jobs)
## 20190412	qsub jobs 7430787-7430792, 7430794-7430797 (10 jobs)
## 20190412 Change EOL conversion in Notepad++ to Unix from default Windows CR LF. Notepad++ > Edit > EOL conversion > Unix. Had errors due to backslash r yesterday.
## 20190411	qsub jobs 7347160-7347162, 7347164-7347174, 7347176-7347186, 7347188-7347198, 7347200-7347208 (discontinuous job numbers)
##---------------------------------------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806"
filePath_jobScript_continuous="${locScripts}/MR_step08-02-04_jobScript_calculate-LDSC-SNP-heritability-for-a-continuous-trait.sh"
filePath_jobScript_binary="${locScripts}/MR_step08-02-05_jobScript_calculate-LDSC-SNP-heritability-for-a-binary-trait.sh"
locHistory="${homeDir}/history"

#-------------------------------------------
# Folders under lunC working
#-------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data"
loc_LDSC="${workingDir}/MR_ICC_GSCAN_201806/LD-score-correlation"
loc_LDSC_input="${loc_LDSC}/input"
loc_LDSC_SNP_heritability="${loc_LDSC}/output/SNP-heritability"
loc_LDSC_SNP_heritability_archive=${loc_LDSC_SNP_heritability}/archive
pbs_output_dir="${loc_LDSC_SNP_heritability}/pbs_output"
pbs_output_dir_archive=${pbs_output_dir}/archive
mkdir -p ${loc_LDSC_SNP_heritability} ${pbs_output_dir} ${loc_LDSC_SNP_heritability_archive} ${pbs_output_dir_archive}

# Location of a tsv file to loop thru
filePath_munged_GWAS=${loc_LDSC_input}/file-info_munged-QCed-GWASs.tsv

#-----------------------------------------------------------------------------------------------------
# If running a reanalysis, first archive existing result files and PBS log files in the output folders 
## Comment out the code after done
#-----------------------------------------------------------------------------------------------------
#find ${loc_LDSC_SNP_heritability} -name "SNP-heritability_*.log" -exec mv {} ${loc_LDSC_SNP_heritability_archive} \;
#find ${pbs_output_dir} -name "h2_*.pbs.*" -exec mv {} ${pbs_output_dir_archive} \;

#-----------------------------------------------------------------------------------------------------
# Get number of lines including header row
#-----------------------------------------------------------------------------------------------------
num_lines=$(awk -F"\t" 'END {print NR}' $filePath_munged_GWAS) # 12

#-----------------------------------------------------------------------------------------------------
# Set up resources requested for submitting PBS jobs
#-----------------------------------------------------------------------------------------------------
num_cpu=1;
runTime_requested=00:10:00; 
memory_requested=5gb;
cpu_brand=Intel; # This avoids the use of node 34, a very slow one from AMD

#-------------------------------------------------------------------------------------------------------------
# Loop through line 2 to last line of a tsv file where each line contains information about a munged GWAS file
## Number of iterations: 11
## Skip first line (header) of the tsv file
#-------------------------------------------------------------------------------------------------------------
count=0
for ((line=2;line <=${num_lines};line++));do
	count=$(($count+1));
	jobName="h2_$count"
	echo "==================================== Iteration $count ======================================================="
	field_1=$(awk '(NR=='$line') {print $1}' $filePath_munged_GWAS); # Path of munged file
	field_2=$(awk '(NR=='$line') {print $2}' $filePath_munged_GWAS); # consortium
	field_3=$(awk '(NR=='$line') {print $3}' $filePath_munged_GWAS); # trait
	field_4=$(awk '(NR=='$line') {print $4}' $filePath_munged_GWAS); # substance type of the trait
	field_5=$(awk '(NR=='$line') {print $5}' $filePath_munged_GWAS); # variable type of the trait
	field_6=$(awk '(NR=='$line') {print $6}' $filePath_munged_GWAS); # Sample size of a continuous trait
	field_9=$(awk '(NR=='$line') {print $9}' $filePath_munged_GWAS); # population prevalence for binary traits
	field_10=$(awk '(NR=='$line') {print $10}' $filePath_munged_GWAS); # sample prevalence for binary traits
	echo "field_1=$field_1";
	echo "field_2=$field_2";
	echo "field_3=$field_3";
	echo "field_4=$field_4";
	echo "field_5=$field_5";
	echo "field_6=$field_6";
	echo "field_9=$field_9";
	echo "field_10=$field_10";
	# Conditionally run the script to calculate SNP heritability
	if [ "$field_5" = "binary" ]; then
		echo "qsub -N $jobName -v v_file_path_munged=${field_1},v_consortium=${field_2},v_trait=${field_3},v_population_prevalence=${field_9},v_sample_prevalence=${field_10},v_folderPath_output=${loc_LDSC_SNP_heritability} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${cpu_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${filePath_jobScript_binary}"
		qsub -N $jobName -v v_file_path_munged=${field_1},v_consortium=${field_2},v_trait=${field_3},v_population_prevalence=${field_9},v_sample_prevalence=${field_10},v_folderPath_output=${loc_LDSC_SNP_heritability} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${cpu_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${filePath_jobScript_binary}
	elif [ "$field_5" = "continuous" ]; then
		echo "qsub -N $jobName -v v_file_path_munged=${field_1},v_consortium=${field_2},v_trait=${field_3},v_folderPath_output=${loc_LDSC_SNP_heritability} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${cpu_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${filePath_jobScript_continuous}"
		qsub -N $jobName -v v_file_path_munged=${field_1},v_consortium=${field_2},v_trait=${field_3},v_folderPath_output=${loc_LDSC_SNP_heritability} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${cpu_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${filePath_jobScript_continuous};
	fi

done

#cp -n ${locScripts}/MR_step08-02-06_jobSubmit_calculate-LDSC-SNP-heritability.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#



