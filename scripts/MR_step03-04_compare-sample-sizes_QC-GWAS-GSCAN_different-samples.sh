#!/bin/bash
## file name: MR_step03-01_QC_GWAS-GSCAN-ICC.sh
## old file name: 
## modified from: PRS_UKB_201711_step01-04_jobSubmi_QC-GWAS.sh
## date created: 20180719
## purpose: Quality Control GWAS files by removing duplicated SNPs, ambiguous SNPs
## Run dependency: 
## How to run this script: line by line

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20181119	QC GSCAN GWAS files (sample excluded ICC, UKB, 23andme, QIMR BLTS)
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input	${loc_no23andMe_QC3}/ai_no23andme.ambiguousSNPRemoved
# Input	${loc_no23andMe_QC3}/cpd_no23andme.ambiguousSNPRemoved
# Input	${loc_no23andMe_QC3}/dpw_no23andme.ambiguousSNPRemoved
# Input	${loc_no23andMe_QC3}/sc_no23andme.ambiguousSNPRemoved
# Input	${loc_no23andMe_QC3}/si_no23andme.ambiguousSNPRemoved

# Input	${loc_noICC_QC3}/ai_noICC.ambiguousSNPRemoved
# Input	${loc_noICC_QC3}/cpd_noICC.ambiguousSNPRemoved
# Input	${loc_noICC_QC3}/dpw_noICC.ambiguousSNPRemoved
# Input	${loc_noICC_QC3}/sc_noICC.ambiguousSNPRemoved
# Input	${loc_noICC_QC3}/si_noICC.ambiguousSNPRemoved

# Outpu ${output}/ai_header_added
# Outpu ${output}/cpd_header_added
# Outpu ${output}/dpw_header_added
# Outpu ${output}/sc_header_added
# Outpu ${output}/si_header_added
#---------------------------------------------------------------------------------------------------------------
## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";
locGSCAN="${homeDir}/LabData/Lab_NickM/lunC/GSCAN";
locICC="${homeDir}/LabData/Lab_NickM/lunC/international_cannabis_consortium";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";
loc_no23andMe="$locMR/no23andMe_results";
loc_noICC="$locMR/noICC_results";

loc_no23andMe_QC3=$loc_no23andMe/QC3_remove_ambiguousSNPs_indel;
loc_noICC_QC3=$loc_noICC/QC3_remove_ambiguousSNPs_indel;

locFileJoiner="/reference/genepi/GWAS_release/Release8/Scripts/FileJoiner"
loc_destination="/mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/collab/lunC/MR_ICC_GSCAN_201806/data/SNPs.compare.no23andMe.noICC"
input="${loc_destination}/input"
output="${loc_destination}/output"
mkdir -p $input $output

# do a left join
#file1=${loc_no23andMe_QC3}/ai_no23andme.ambiguousSNPRemoved
#file2=${loc_noICC_QC3}/ai_noICC.ambiguousSNPRemoved

# Left join no23andMe (left table) and noICC (right table) GWAS from the same trait
## Create a two column file where column1 is file path of no23andMe and column2 is file path of noICC
paste $loc_no23andMe_QC3/filePath $loc_noICC_QC3/filePath > ${input}/filePath

IFS=$'\n'; 
count=0;
for line in $(cat ${input}/filePath); do
	count=$(($count+1))
	filePath_left_table=$(echo ${line} | awk '{print $1}')
	filePath_righ_table=$(echo ${line} | awk '{print $2}')
	traitName_left=$(basename $filePath_left_table | cut -d"_" -f1 )
	traitName_righ=$(basename $filePath_righ_table | cut -d"_" -f1 )
	echo "=============================================== iteration $count ================================="
	echo "filePath_left_table=$filePath_left_table"
	echo "filePath_righ_table=$filePath_righ_table"
	echo "traitName_left=$traitName_left"
	if [ "$traitName_left" = "$traitName_righ" ]; then
		# Perform an inner join
		$locFileJoiner -quiet $filePath_left_table,headerline,sep=tabs,key=3,outfields=3,10-12 $filePath_righ_table,headerline,sep=tabs,key=3,outfields=10-12 > ${output}/${traitName_left}
		# Add a header to the joined file that is tab separated
		cat <(echo -e "RSID\tN_no23andMe\tEFFECTIVE_N_no23andMe\tNumber_of_Studies_no23andMe\tN_noICC\tEFFECTIVE_N_noICC\tNumber_of_Studies_noICC") <(tail -n +2  ${output}/${traitName_left}) > ${output}/${traitName_left}_header_added
	else
		echo "Non-matched trait! No file joining performed"
	fi
done

##---------------------------------This is the end of this file-------------------------------------------------##