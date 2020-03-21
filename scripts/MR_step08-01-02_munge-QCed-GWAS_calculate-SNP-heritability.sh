#!/bin/bash
 
## File name: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step08-01-02_munge-QCed-GWAS_calculate-SNP-heritability.sh
## Modified from: MR_step08-01_calculate-genetic-correlations-between-exposure-outcome-GWAS_LD-score-regression.sh PRS_UKB_201711_step22-01_calculate-genetic-correlations-between-discovery-phenotypes_LD-score-regression.sh
## Date created: 20190410
## Note: 
## Purpose: Munge QCed GWAS files and calculate SNP heritability using LD score correlation python scripts
## Run dependency: 
### (1) ${locScripts}/MR_step08-01-01_create-files-to-loop-through-running-LDSC-munge.R
## How to run this file: 
### . $locScripts/MR_step08-01-02_munge-QCed-GWAS_calculate-SNP-heritability.sh > $locHistory/run-shell-script_20190813_load-module-ldsc_munge-QCed-GWAS-files
## Note: Don't open the CSV files by right clicking edit with Excel. This may cause all columns appear in single column. Instead, import the CSV files in Excel locally, edit them and save them again as CSV, and transfer them back to remote folder. See instruction:  https://support.goteamup.com/hc/en-us/articles/203460681-FAQ-When-opening-CSV-file-all-data-appears-in-one-column

# Type		Files 
#------------------------------------------------------------------------------------------------------
# Input 	/reference/data/UKBB_500k/versions/lab_stuartma/LDSC_res/1000HGP/w_hm3.snplist	
# Input 	module load ldsc munge_sumstats.py
# Input		${loc_LDSC_input}/file-info_QCed-GWASs.tsv

# Outpu		${loc_LDSC_input}/HapMap3-SNPs-found-in_QCed-GWAS-* (11 files)
# Outpu		${loc_LDSC_input}/header-fieldNumber-colname_* (11 files)
# Outpu		${loc_LDSC_output}/HapMap3-SNPs-found-in_QCed-GWAS-*_munged.sumstats.gz (11 files)
# Outpu		${loc_LDSC_output}/HapMap3-SNPs-found-in_QCed-GWAS-*_munged.log (11 files)
#------------------------------------------------------------------------------------------------------
## Time 	Change
##---------------------------------------------------------------------------------------------------------------------
## 20190813	Munged 11 GWAS files above (used module load ldsc rather than ${ldsc_dir}/munge_sumstats.py)
## 20190411	Munged 10 GWAS files above
##---------------------------------------------------------------------------------------------------------------------

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

#-------------------------------------------
# Folders under lunC working
#-------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR=${workingDir}/MR_ICC_GSCAN_201806/data
loc_LDSC=${workingDir}/MR_ICC_GSCAN_201806/LD-score-correlation
loc_LDSC_input=${loc_LDSC}/input
loc_LDSC_output=${loc_LDSC}/output/munged-GWASs
loc_LDSC_output_archive=${loc_LDSC_output}/archive
mkdir -p ${loc_LDSC_output} ${loc_LDSC_output_archive}

# Location of a tsv file to loop thru
filePath_GWAS_info=${loc_LDSC_input}/file-info_QCed-GWASs.tsv

#-------------------------------------------
# Archive existing output files before a reanalysis
#-------------------------------------------
# find ${loc_LDSC_output} -name "HapMap3-SNPs-found-in_QCed-GWAS-*" -exec mv {} ${loc_LDSC_output_archive} \;

#-------------------------------------------------------------------
# Loop thru a tsv file skipping its first row (header)
#-------------------------------------------------------------------
IFS=$'\n';
count=0;
for line in $(tail -n+2 $filePath_GWAS_info); do
#for line in $(sed -n '8p' $filePath_GWAS_info); do
	# In every line, get values of every field
	filePath=$(echo $line | cut -d$'\t' -f1)
	colname_SNP=$(echo $line | cut -d$'\t' -f2)
	colname_effect_allele=$(echo $line | cut -d$'\t' -f3)
	colname_other_allele=$(echo $line | cut -d$'\t' -f4)
	colname_beta=$(echo $line | cut -d$'\t' -f5)
	colname_SE=$(echo $line | cut -d$'\t' -f6)
	colname_p_value=$(echo $line | cut -d$'\t' -f7)
	consortium=$(echo $line | cut -d$'\t' -f8)
	trait=$(echo $line | cut -d$'\t' -f9)
	trait_type=$(echo $line | cut -d$'\t' -f10)
	N_continuous=$(echo $line | cut -d$'\t' -f11)
	N_cases=$(echo $line | cut -d$'\t' -f12)
	N_controls=$(echo $line | cut -d$'\t' -f13)
	population_prevalence=$(echo $line | cut -d$'\t' -f14)
	file_delimiter=$(echo $line | cut -d$'\t' -f15)
	if [ "$file_delimiter" = "tab" ]; then 
	file_delimiter_symbol="\t";
	elif [ "$file_delimiter" = "space" ]; then 
	file_delimiter_symbol=" ";
	fi
	#--------------------------------------------------------------------
	# Export field number and column name of input file's header row as a 2 column file
	#--------------------------------------------------------------------
	head -1 $filePath | awk '{for (i=1; i <=NF; i++) print i, $i}' > ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait};
	# Return the field position of individual column names
	position_SNP=$(awk -v pattern=$colname_SNP '{if ($2 ~ "^" pattern "$") print $1}' ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait});
	position_effect_allele=$(awk -v pattern=$colname_effect_allele '{if ($2 ~ "^" pattern "$") print $1}' ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait})
	position_other_allele=$(awk -v pattern=$colname_other_allele '{if ($2 ~ "^" pattern "$") print $1}' ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait})
	position_beta=$(awk -v pattern=$colname_beta '{if ($2 ~ "^" pattern "$") print $1}' ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait})
	position_SE=$(awk -v pattern=$colname_SE '{if ($2 ~ "^" pattern "$") print $1}' ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait})
	position_p_value=`awk -v pattern=${colname_p_value} '{if ($2 ~ "^" pattern "$") print $1}' ${loc_LDSC_input}/header-fieldNumber-colname_${consortium}-${trait}`
	#--------------------------------------------------------------------
	# Inspect assigned variables to see if they are correct
	#--------------------------------------------------------------------
	echo "filePath=$filePath"; 
	echo "colname_SNP=$colname_SNP, colname number=$position_SNP";
	echo "colname_effect_allele=$colname_effect_allele, colname number=$position_effect_allele";
	echo "colname_other_allele=$colname_other_allele, colname number=$position_other_allele";
	echo "colname_beta=$colname_beta, colname number=$position_beta";
	echo "colname_SE=$colname_SE, colname number=$position_SE";
	echo "colname_p_value=$colname_p_value, colname number=$position_p_value";
	echo "consortium=$consortium"; 
	echo "trait=$trait"; 
	echo "trait_type=$trait_type";
	echo "N_continuous=$N_continuous";
	echo "N_cases=$N_cases";
	echo "N_controls=$N_controls";
	echo "population_prevalence=$population_prevalence";
	#---------------------------------------------------------------------------------------------------------------------------------------------
	# Find reference SNP lists in the input file, exporting only 6 columns using awk
	## {a[$1];next;} Compare file2 ${filePath} with field 1 of file1 ${ref_snplist}
	## {if($field_1 in a ) print $field_1,$field_2,$field_3,$field_4,$field_5,$field_6} If SNP column of file 2 is found in file1, then output the 6 selected columns from file 2
	#---------------------------------------------------------------------------------------------------------------------------------------------
	 awk -F "${file_delimiter_symbol}" -v field_1=${position_SNP} -v field_2=${position_effect_allele} -v field_3=${position_other_allele} -v field_4=${position_beta} -v field_5=${position_SE} -v field_6=${position_p_value} 'BEGIN {print "SNP","effect_allele","other_allele","beta","SE","p_value"}(NR==FNR)&&(NR!=1) {a[$1];next;} {if($field_1 in a ) print $field_1,$field_2,$field_3,$field_4,$field_5,$field_6}' ${ref_snplist} ${filePath} > ${loc_LDSC_input}/HapMap3-SNPs-found-in_QCed-GWAS-${consortium}-${trait};
	#--------------------------------------------------------------------------------------------------------------------------------------------
	# Conditionally munge GWAS files in LDSC python script
	## if the trait of GWAS is continuous, run the first then action
	## if the trait of GWAS is binary, run the other then action
	#---------------------------------------------------------------------------------------------------------------------------------------------
	module load ldsc
	if [ "$trait_type" = "continuous" ]; 
	then munge_sumstats.py --sumstats ${loc_LDSC_input}/HapMap3-SNPs-found-in_QCed-GWAS-${consortium}-${trait} --N ${N_continuous} --out ${loc_LDSC_output}/HapMap3-SNPs-found-in_QCed-GWAS-${consortium}-${trait}_munged --merge-alleles ${ref_snplist};
	elif [ "$trait_type" = "binary" ]; 
	then munge_sumstats.py --sumstats ${loc_LDSC_input}/HapMap3-SNPs-found-in_QCed-GWAS-${consortium}-${trait} --N-cas ${N_cases} --N-con ${N_controls} --out ${loc_LDSC_output}/HapMap3-SNPs-found-in_QCed-GWAS-${consortium}-${trait}_munged --merge-alleles ${ref_snplist};
	fi
done

#cp -n ${locScripts}/MR_step08-01-02_munge-QCed-GWAS_calculate-SNP-heritability.sh ${locScripts}/MR_step08-02-02_calculate-genetic-correlation-LDSC.sh
#-------------------------------------------------------------------------------------------------#
#--------------- This is the end of this file ----------------------------------------------------#
#-------------------------------------------------------------------------------------------------#



