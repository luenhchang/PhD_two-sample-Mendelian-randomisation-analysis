#!/bin/bash
 
## File name: ${locScripts}/MR_step08-02-03_jobSubmit_calculate-genetic-correlation-LDSC.sh
## Modified from: ${locScripts}/MR_step08-01-02_munge-QCed-GWAS_calculate-SNP-heritability.sh 
## Date created: 20190411
## Note: 
## Purpose: calculate genetic correlation between two munged QCed GWAS files that are different traits and consortium using LD score correlation python scripts
## Run dependency: 
### (1) ${locScripts}/MR_step08-01-01_create-files-to-loop-through-running-LDSC-munge.R
### (2) ${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh
### (3) module load ldsc munge_sumstats.py

## How to run this file: 
### . ${locScripts}/MR_step08-02-03_jobSubmit_calculate-genetic-correlation-LDSC.sh > ${locHistory}/run-shell-script_20190813_calculate-rG-between-munge-QCed-GWAS-files
## Note: Don't open the CSV files by right clicking edit with Excel. This may cause all columns appear in single column. Instead, import the CSV files in Excel locally, edit them and save them again as CSV, and transfer them back to remote folder. See instruction:  https://support.goteamup.com/hc/en-us/articles/203460681-FAQ-When-opening-CSV-file-all-data-appears-in-one-column

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input 	/reference/data/UKBB_500k/versions/lab_stuartma/LDSC_res/1000HGP/w_hm3.snplist	
# Input		${loc_LDSC_input}/file-info_munged-QCed-GWASs.tsv

# Outpu		${loc_LDSC_rG}/genetic-correlation_between_*-*_and_*-*.log (45)
#------------------------------------------------------------------------------------------------------
## Time 	Change
##---------------------------------------------------------------------------------------------------------------------
## 20190813	qsub jobs 8616004-8616059 (55 jobs)
## 20190412 Change EOL conversion in Notepad++ to Unix from default Windows CR LF. Notepad++ > Edit > EOL conversion > Unix. Had errors due to backslash r yesterday.
## 20190411	qsub jobs 7347160-7347162, 7347164-7347174, 7347176-7347186, 7347188-7347198, 7347200-7347208 (discontinuous job numbers)
##---------------------------------------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806"
filePath_jobScript="${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh"
locHistory="${homeDir}/history"

#----------------------------------------------
# Folder locations under lunC working directory
#----------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data"
loc_LDSC="${workingDir}/MR_ICC_GSCAN_201806/LD-score-correlation"
loc_LDSC_input="${loc_LDSC}/input"
loc_LDSC_munged="${loc_LDSC}/output/munged-GWASs"
loc_LDSC_rG="${loc_LDSC}/output/genetic-correlations"
loc_LDSC_rG_archive=${loc_LDSC_rG}/archive
pbs_output_dir="${loc_LDSC_rG}/pbs_output"
pbs_output_dir_archive=${pbs_output_dir}/archive
mkdir -p ${loc_LDSC_rG} ${pbs_output_dir} ${loc_LDSC_rG_archive} ${pbs_output_dir_archive}

# Location of a tsv file to loop thru
filePath_munged_GWAS=${loc_LDSC_input}/file-info_munged-QCed-GWASs.tsv

#-----------------------------------------------------------------------------------------------------
# If running a reanalysis, first archive existing result files and PBS log files in the output folders 
## Comment out the code after done
#-----------------------------------------------------------------------------------------------------
#find ${loc_LDSC_rG} -name "genetic-correlation_between_*-*_and_*-*.log" -exec mv {} ${loc_LDSC_rG_archive} \;
#find ${pbs_output_dir} -name "rG_iteration_*" -exec mv {} ${pbs_output_dir_archive} \;

#----------------------------------------------
# Get number of lines including header row
#----------------------------------------------
num_lines=$(awk -F"\t" 'END {print NR}' $filePath_munged_GWAS) # 11

#----------------------------------------------
# Set up resources requested for submitting PBS jobs
#----------------------------------------------
num_cpu=1;
runTime_requested=3:00:00; 
memory_requested=2gb;
cpu_brand=Intel; # This avoids the use of node 34, a very slow one from AMD

#-------------------------------------------------------------------
# Loop through unique combinations of 11 munged files
## Number of iterations: 55 (11 choose 2) 
#-------------------------------------------------------------------
# Iterate thru 55 combinations 
## Skip first line (header)
count=0
for ((i=2;i <=${num_lines}-1;i++));do
	for ((j=$i+1;j<=${num_lines};j++));do
		count=$(($count+1));
		jobName="rG_iteration_$count"
		echo "==================================== Iteration $count ======================================================="
		#---------------------------------------------------------------------------------------------------
		# Assign values of field 1-3 per two line counter variables $i and $j of the input file to variables
		#---------------------------------------------------------------------------------------------------
		file=${filePath_munged_GWAS} 
		mungedFile_1_field_1=$(awk '(NR=='$i') {print $1}' $file); # file path
		mungedFile_1_field_2=$(awk '(NR=='$i') {print $2}' $file); # consortium
		mungedFile_1_field_3=$(awk '(NR=='$i') {print $3}' $file); # trait
		mungedFile_1_field_9=$(awk '(NR=='$i') {print $9}' $file); # population prevalence for binary traits
		mungedFile_1_field_10=$(awk '(NR=='$i') {print $10}' $file); # sample prevalence for binary traits

		mungedFile_2_field_1=$(awk '(NR=='$j') {print $1}' $file); # file path
		mungedFile_2_field_2=$(awk '(NR=='$j') {print $2}' $file); # consortium
		mungedFile_2_field_3=$(awk '(NR=='$j') {print $3}' $file); # trait
		mungedFile_2_field_9=$(awk '(NR=='$j') {print $9}' $file); # population prevalence for binary traits
		mungedFile_2_field_10=$(awk '(NR=='$j') {print $10}' $file); # sample prevalence for binary traits

		echo "munged file 1: ${mungedFile_1_field_2}-${mungedFile_1_field_3}, ${mungedFile_1_field_1}, ${mungedFile_1_field_9}, ${mungedFile_1_field_10}"
		echo "munged file 2: ${mungedFile_2_field_2}-${mungedFile_2_field_3}, ${mungedFile_2_field_1}, ${mungedFile_2_field_9}, ${mungedFile_2_field_10}"
		#-------------------------------------------------------------------------------------------------------------
		# Calculate genetic correlation between munged file 1 and munged file 2 where the 2 files are different traits
		#-------------------------------------------------------------------------------------------------------------
		echo "qsub -N $jobName -v v_MF_1_filePath=${mungedFile_1_field_1},v_MF_1_consortium=${mungedFile_1_field_2},v_MF_1_trait=${mungedFile_1_field_3},v_MF_1_population_prevalence=${mungedFile_1_field_9},v_MF_1_sample_prevalence=${mungedFile_1_field_10},v_MF_2_filePath=${mungedFile_2_field_1},v_MF_2_consortium=${mungedFile_2_field_2},v_MF_2_trait=${mungedFile_2_field_3},v_MF_2_population_prevalence=${mungedFile_2_field_9},v_MF_2_sample_prevalence=${mungedFile_2_field_10},v_folderPath_output=${loc_LDSC_rG} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${cpu_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${filePath_jobScript}"
		qsub -N $jobName -v v_MF_1_filePath=${mungedFile_1_field_1},v_MF_1_consortium=${mungedFile_1_field_2},v_MF_1_trait=${mungedFile_1_field_3},v_MF_1_population_prevalence=${mungedFile_1_field_9},v_MF_1_sample_prevalence=${mungedFile_1_field_10},v_MF_2_filePath=${mungedFile_2_field_1},v_MF_2_consortium=${mungedFile_2_field_2},v_MF_2_trait=${mungedFile_2_field_3},v_MF_2_population_prevalence=${mungedFile_2_field_9},v_MF_2_sample_prevalence=${mungedFile_2_field_10},v_folderPath_output=${loc_LDSC_rG} -l ncpus=${num_cpu},walltime=${runTime_requested},mem=${memory_requested},chip=${cpu_brand} -e ${pbs_output_dir}/${jobName}.pbs.err -o ${pbs_output_dir}/${jobName}.pbs.out ${filePath_jobScript}; 
	done
done

#cp -n ${locScripts}/MR_step08-02-03_jobSubmit_calculate-genetic-correlation-LDSC.sh ${locScripts}/MR_step08-02-06_jobSubmit_calculate-LDSC-SNP-heritability.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#



