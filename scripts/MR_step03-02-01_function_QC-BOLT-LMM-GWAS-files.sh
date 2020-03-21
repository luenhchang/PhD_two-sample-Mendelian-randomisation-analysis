#!/bin/bash
## File path: /mnt/backedup/home/lunC/scripts/Bash_functions/quality-control-GWAS-files.sh 
## Modified from: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step03-02_QC_GWAS-UKB.sh
## Date created: 20190811
## Purpose: Quality Control GWAS files by removing duplicated SNPs, ambiguous SNPs
## Run dependency: 
## Functions internal: copy_files_get_filePaths
## References: [bash: how to pass command line arguments containing special characters](https://superuser.com/questions/163515/bash-how-to-pass-command-line-arguments-containing-special-characters)
## [How to call bash functions](https://superuser.com/questions/106272/how-to-call-bash-functions)
## How to run this script: 

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input 
# Outpu 
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------
# Folder locations under my home
#---------------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locBashFunctions="${homeDir}/scripts/Bash_functions"
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";

#------------------------------------------------------------------
# Create functions to test function calls
#------------------------------------------------------------------
# A function without arguments
function test1 {
	echo "testing function test1";
}

# An example of calling the function above
# source /mnt/backedup/home/lunC/scripts/Bash_functions/quality-control-GWAS-files.sh; 
# test1;

# A function with 1 argument
function test2 {
	local var1=${1}
	echo "testing function test2, passed argument is $var1 ";
}

# source /mnt/backedup/home/lunC/scripts/Bash_functions/quality-control-GWAS-files.sh; 
# test2 "my argument 1";

#------------------------------------------------------------------
# Create a function to copy GWAS files and get file paths
#------------------------------------------------------------------
function copy_files_get_filePaths {
	# Pass arguments to local shell variables within this function. Default is a global variable in shell. 
	local source_folder_path=$1;
	local source_file_names=$2;	
	local destin_folder_path=$3;
	# Print the shell variable values
	echo "source_folder_path=${source_folder_path}";
	echo "source_file_name=${source_file_name}";
	echo "destin_folder_path=${destin_folder_path}";
	# Execute commands using the shell variables
	## Copy GWAS files from source folder to destination folder	
	find ${source_folder_path} -name "${source_file_name}" -exec cp {} ${destin_folder_path} \; 
	
	# Store file paths in a file 
	#echo "realpath ${destin_folder_path}/${source_file_name} > ${destin_folder_path}/filePath"
	realpath ${destin_folder_path}/${source_file_name} > ${destin_folder_path}/filePath;
}

# A working example of calling the function 
# copy_files_get_filePaths "${ref_UKB_caffeine_NA20453_GWAS}" "revised_bolt_imputed_ukb_imp_chr*v3_caffeine.per.day.bgen.assoc" "$locUKB_caffeine_raw"

#-------------------------------------------------------------------------------------------------------------
# Create a function to loop thru every line of an input file. For each line (a tab-separated file path ), print lines with more than 1 occurrence of SNP rs number using awk two file processing
#-------------------------------------------------------------------------------------------------------------
function remove_multiple_occurrences_SNPs {
	local field_RSNum=${1}; # Filed position of SNP Rs ID
	local input_file_path=${2};
	local output_folder_path_step1=${3}; # Folder that holds duplicate lines of GWAS files
	local output_folder_path_step2=${4}; # Folder where GWAS with only 1 occurrence of SNP ID are exported
	echo "field_RSNum=${field_RSNum}"
	echo "input_file_path=${input_file_path}"
	echo "output_folder_path_step1=${output_folder_path_step1}"
	
	# Loop thru every line of an input file, where each line is a file path  
	rm -f ${output_folder_path_step2}/filePath; # Force remove the output file if it is already there
	for filePath in `cat $input_file_path`; do
		local fileName=`basename $filePath`;
		echo "fileName=$fileName";
		# Print lines with more than 1 occurrence of SNP rs number using awk two file processing 
		echo "awk -F\"\t\" -v fieldPosSNPID=\$field_RSNum 'NR==FNR {count[\$fieldPosSNPID]++;next} count[\$fieldPosSNPID]>1' ${filePath} ${filePath} > ${output_folder_path_step1}/${fileName}.allOccurrenceDup"
		awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' ${filePath} ${filePath} > ${output_folder_path_step1}/${fileName}.allOccurrenceDup ;
		
		# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
		## output files are headerless		
		awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > ${output_folder_path_step2}/${fileName}.allOccurrenceDupRemoved;
		# Store paths of GWAS files with unique SNPs in a file
		realpath ${output_folder_path_step2}/${fileName}.allOccurrenceDupRemoved >> ${output_folder_path_step2}/filePath
	done;
}

# A working example of calling the function
# remove_multiple_occurrences_SNPs 1 "$locUKB_caffeine_raw/filePath" "$locUKB_caffeine_QC1" "$locUKB_caffeine_QC2"

#-------------------------------------------------------------------------------------------------------------
# Create a function to loop through an input file, where each line is a tab-separated file path. For each file, 
## (1) further remove non SNPid (rs number) in the SNP field from files from previous function
## (2) remove ambiguous SNPs and insertions, deletions
## The input files have 16 columns. Make sure headers in the awk BEGIN match the number
#-------------------------------------------------------------------------------------------------------------
function remove_nonSNPid_ambiguous_SNPs_insertions_deletions {
	local input_file_path=${1};
	local output_folder_path=${2};
	local output_file_name=${3};
	echo "input_file_path=${input_file_path}"
	echo "output_folder_path=${output_folder_path}"
	echo "output_file_name=${output_file_name}"
	# Remove lines with non SNPid (rs number), ambiguous SNPs and insertions, deletions
	rm -f ${output_folder_path}/filePath # Force remove the output file if it is already there
	rm -f ${output_folder_path}/new-header
	for filePath in `cat $input_file_path`;do 
		local fileName=`basename $filePath`; 
		echo "fileName=$fileName";
		awk 'BEGIN {FS="\t";OFS="\t"; print "SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM",
		"P_BOLT_LMM" } {if( $1 ~/rs/) print $0}' ${filePath} | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${output_folder_path}/${fileName}.ambiguousSNPRemoved;
		# Save file paths in a file
		realpath ${output_folder_path}/${fileName}.ambiguousSNPRemoved >> ${output_folder_path}/filePath
		# Make new headers in a file. Change P_BOLT_LMM_INF to PVALUE for LD score regression to interpret the p value column easily
		head -1 ${output_folder_path}/${fileName}.ambiguousSNPRemoved | sed 's/P_BOLT_LMM/PVALUE/g' >> ${output_folder_path}/new-header
	done;
	
	# Vertically combine all 22 files excluding their header rows
	cat ${output_folder_path}/filePath | xargs awk '(FNR>1){print $0}'> ${output_folder_path}/${output_file_name}
	# Write header back to the file
	cd ${output_folder_path};
	cat <(head -1 new-header) <(cat ${output_file_name}) > ${output_file_name}_headed

}

# A working example of calling the function
#remove_nonSNPid_ambiguous_SNPs_insertions_deletions "$locUKB_caffeine_QC2/filePath" "$locUKB_caffeine_QC3" "QCed-GWAS-UKB-caffeine-consumed-per-day"

#cp -n /mnt/backedup/home/lunC/scripts/Bash_functions/quality-control-GWAS-files.sh /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step03-02-01_function_QC-BOLT-LMM-GWAS-files.sh
