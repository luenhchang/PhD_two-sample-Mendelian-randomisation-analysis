#!/bin/bash

## File path: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step05-01_merge_LD-clumped-SNP_GWAS.sh
## old file name: zMR_step05_merge_LD-clumped-SNP_GWAS.sh
## modified from: 
## date created: 20190404
## purpose: leftJoin SNPs clumped at MR_step04 (left table) and GWAS files QCed at MR_step03 (right table)
## Run dependency:  /reference/genepi/GWAS_release/Release8/Scripts/FileJoiner
## How to run this script: . ${locScripts}/MR_step05_merge_LD-clumped-SNP_GWAS.sh > ${locHistory}/jobSubmitted_20190404_leftJoin_LD-clumped-SNPs_and_QCed-GWAS-files-GSCAN-UKB-ICC

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20190812	LeftJoin clumped GWASs and QCed GWASs on UKB caffeine consumed per day in people without data in 20453.
## 20190408 LeftJoin clumped GWASs and QCed GWASs from GSCAN. Created 10 output files. Previous run didn't deal with inconsistent header (13 fields)and data columns (14 fields)
## 20190404 LeftJoin clumped GWASs and QCed GWASs. Created 20 output files
## 20190306 LeftJoin clumped GSCAN GWAS (LD window=10000kb) and QCed GSCAN GWAS. Created 10 output files
## 20190305 LeftJoin clumped GSCAN GWAS (LD window=500kb) and QCed GSCAN GWAS. Created 10 output files
## 20181122	left join files for GSCAN. Created 10 output files
## 20180823 left join files for UKB ESDPW, CCPD, PYOS in people without data in 20453. Created files 4*3
## 20180822 left join files for UKB ESDPW in people without data in 20453. Created 4 output files
## 20180806	left join files for UKB 3456 in people without data in 20345. Created 4 output files
## 20180729 left join files for GSCAN. Created 20 output files
## 20180719	left join files for GSCAN
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# [Clumped SNPs]
# Input ${locGSCAN_LDOut}/*_noICC_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (5 traits * 2 Clumping criteria= 10 files)
# Input ${locUKB3456_LDOut}/GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Input $locUKB_ESDPW_LDOut/GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files) 
# Input $locUKB_CCPD_LDOut/GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Input $locUKB20161_LDOut/GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Input $locICC_LDOut/GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Input $locUKB_caffeine_LDOut/GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)

# [QCed GWAS files]
# Input ${locGSCAN_QC3}/ai_noICC.ambiguousSNPRemoved
# Input ${locGSCAN_QC3}/cpd_noICC.ambiguousSNPRemoved
# Input ${locGSCAN_QC3}/dpw_noICC.ambiguousSNPRemoved
# Input ${locGSCAN_QC3}/sc_noICC.ambiguousSNPRemoved
# Input ${locGSCAN_QC3}/si_noICC.ambiguousSNPRemoved
# Input ${locUKB3456_QC3}/QCed-GWAS-UKB3456_headed
# Input $locUKB_ESDPW_QC3/QCed-GWAS-UKB-ESDPW_headed
# Input $locUKB_CCPD_QC3/QCed-GWAS-UKB-CCPD_headed
# Input $locUKB20161_QC3/QCed-GWAS-UKB-PYOS_headed
# Input $locICC/Cannabis_ICC_UKB_small.txt
# Input ${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed

# Outpu realpath $locGSCAN_QC4/GWAS_from-clumped-SNPs_*_noICC_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (10 files)
# Outpu realpath $locUKB3456_QC4/GWAS_from-clumped-SNPs_GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (2 files)
# Outpu realpath $locUKB_ESDPW_QC4/GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (2 files)
# Outpu realpath $locUKB_CCPD_QC4/GWAS_from-clumped-SNPs_GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (2 files)
# Outpu realpath $locUKB20161_QC4/GWAS_from-clumped-SNPs_GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (2 files)
# Outpu realpath $locICC_QC4/GWAS_from-clumped-SNPs_GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (2 files) 
# Outpu realpath $locUKB_caffeine_QC4/GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-*_p2-* (2 files) 
#---------------------------------------------------------------------------------------------------------------
## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
jobScriptFilePath="${locScripts}/MR_step04-01_jobScript_LD-based-SNP-clumping.sh";
locHistory="${homeDir}/history";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";

# Folder locations related to ICC cannabis initiation
locICC="$locMR/ICC-cannabis-ever";
locICC_LDBasedClumping="$locICC/LD-based_SNP_clumping";
locICC_LDOut=$locICC_LDBasedClumping/output;
locICC_QC4=$locICC/QC4_GWAS_from_clumped_SNPs;
mkdir -p $locICC_QC4

# Folders relevant to GSCAN
locGSCAN="$locMR/noICC_results";
locGSCAN_QC3=$locGSCAN/QC3_remove_ambiguousSNPs_indel;
locGSCAN_QC4=$locGSCAN/QC4_GWAS_from_clumped_SNPs;
locGSCAN_LDBasedClumping=$locGSCAN/LD-based_SNP_clumping;
locGSCAN_LDOut=$locGSCAN_LDBasedClumping/output;

# Folder locations related to UKB3456 exluding people with ever using cannabis data
locUKB3456="$locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis";
locUKB3456_QC3=$locUKB3456/QC3_remove_ambiguousSNPs_indel;
locUKB3456_QC4=$locUKB3456/QC4_GWAS_from_clumped_SNPs;
locUKB3456_LDBasedClumping=$locUKB3456/LD-based_SNP_clumping;
locUKB3456_LDOut=$locUKB3456_LDBasedClumping/output;

# Folder locations related to UKB ESDPW exluding people with ever using cannabis data
locUKB_ESDPW=$locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis;
locUKB_ESDPW_QC3=$locUKB_ESDPW/QC3_remove_ambiguousSNPs_indel;
locUKB_ESDPW_QC4=$locUKB_ESDPW/QC4_GWAS_from_clumped_SNPs;
locUKB_ESDPW_LDBasedClumping=$locUKB_ESDPW/LD-based_SNP_clumping;
locUKB_ESDPW_LDOut=$locUKB_ESDPW_LDBasedClumping/output;

# Folder locations related to UKB CCPD exluding people with ever using cannabis data
locUKB_CCPD=$locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis;
locUKB_CCPD_QC3=$locUKB_CCPD/QC3_remove_ambiguousSNPs_indel;
locUKB_CCPD_QC4=$locUKB_CCPD/QC4_GWAS_from_clumped_SNPs;
locUKB_CCPD_LDBasedClumping=$locUKB_CCPD/LD-based_SNP_clumping;
locUKB_CCPD_LDOut=$locUKB_CCPD_LDBasedClumping/output;

# Folder locations related to UKB 20161 PYOS exluding people with ever using cannabis data
locUKB20161=$locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis;
locUKB20161_QC3=$locUKB20161/QC3_remove_ambiguousSNPs_indel;
locUKB20161_QC4=$locUKB20161/QC4_GWAS_from_clumped_SNPs;
locUKB20161_LDBasedClumping=$locUKB20161/LD-based_SNP_clumping;
locUKB20161_LDOut=$locUKB20161_LDBasedClumping/output;

# Folder locations related to UKB estimated caffeine consumed per day
locUKB_caffeine=${locMR}/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis;
locUKB_caffeine_QC3=${locUKB_caffeine}/QC3_remove_ambiguousSNPs_indel;
locUKB_caffeine_QC4=${locUKB_caffeine}/QC4_GWAS_from_clumped_SNPs;
locUKB_caffeine_LDBasedClumping=${locUKB_caffeine}/LD-based_SNP_clumping ;
locUKB_caffeine_LDOut=${locUKB_caffeine_LDBasedClumping}/output;

#mkdir -p $locGSCAN_QC4 $locGSCAN_QC4_LDWindow500 $locGSCAN_QC4_LDWindow10000 $locGSCAN_QC4_LDWindow500_archive $locGSCAN_QC4_LDWindow10000_archive ;
mkdir -p $locGSCAN_QC4 $locGSCAN_QC4_LDWindow500 $locGSCAN_QC4_LDWindow10000 ;
mkdir -p $locUKB3456_QC4 $locUKB_ESDPW_QC4 $locUKB_CCPD_QC4 $locUKB20161_QC4;
mkdir -p ${locUKB_caffeine_QC4}

locFileJoiner="/reference/genepi/GWAS_release/Release8/Scripts/FileJoiner"

#-----------------------------------------------------------------------------------------------------------------
# Archive files
#-----------------------------------------------------------------------------------------------------------------
## Go to the destination folder. Move everything there to the archive folder except for the archive folder itself
#cd $locGSCAN_QC4_LDWindow500;
#mv !(archive_prior20180729) ${locGSCAN_QC4_LDWindow500_archive}/ ;

#cd $locGSCAN_QC4_LDWindow10000;
#mv !(archive_prior20180729) $locGSCAN_QC4_LDWindow10000_archive/ ;

#-----------------------------------------------------------------------------------------------------------------
# Create a file with columns. Utilising the full listing function to extract file paths as one column
## Column 1: LD Clumped SNPs at MR_step04, used as left tables
## Column 2: QCed GWAS files of the same consortia and traits, used as right tables
## Column 3: field number of SNP in QCed GWAS files
### Merging key SNP field number
###-----------------------------
### Consortium 	File1	File2
### GSCAN		$3		$3
### UKB			$3		$1
### ICC			$3		$1
###-----------------------------
## Column 4: Delimiters in QCed GWAS files
## Column 5: Folder paths for left-joined output files
#-----------------------------------------------------------------------------------------------------------------
# List clumped SNP file paths (22 file paths)
realpath ${locGSCAN_LDOut}/*_noICC_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped ${locUKB3456_LDOut}/GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped ${locUKB_ESDPW_LDOut}/GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped ${locUKB_CCPD_LDOut}/GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped ${locUKB20161_LDOut}/GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped ${locICC_LDOut}/GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped ${locUKB_caffeine_LDOut}/GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped > ${locMR}/filePath_clumped-SNPs # wc -l $locMR/filePath_clumped-SNPs 22 files

# List corresponding QCed GWAS file paths (11 file paths)
realpath ${locGSCAN_QC3}/ai_noICC.ambiguousSNPRemoved ${locGSCAN_QC3}/cpd_noICC.ambiguousSNPRemoved ${locGSCAN_QC3}/dpw_noICC.ambiguousSNPRemoved ${locGSCAN_QC3}/sc_noICC.ambiguousSNPRemoved ${locGSCAN_QC3}/si_noICC.ambiguousSNPRemoved ${locUKB3456_QC3}/QCed-GWAS-UKB3456_headed ${locUKB_ESDPW_QC3}/QCed-GWAS-UKB-ESDPW_headed ${locUKB_CCPD_QC3}/QCed-GWAS-UKB-CCPD_headed ${locUKB20161_QC3}/QCed-GWAS-UKB-PYOS_headed ${locICC}/Cannabis_ICC_UKB_small.txt ${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed > ${locMR}/filePath_QCed-GWAS-files # wc -l $locMR/filePath_QCed-GWAS-files 10

# Field number of SNP in QCed GWAS files (10 lines)
(cat <<- _EOF_
3 tabs
3 tabs
3 tabs
3 tabs
3 tabs
1 tabs
1 tabs
1 tabs
1 tabs
1 whitespace
1 tabs
_EOF_
) > ${locMR}/field-number-SNP_delimiters_QCed-GWAS-files # 11 lines


# Write output folder paths to a file (22 lines)
nl=$'\n'
(cat <<- _EOF_
${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}${nl}${locGSCAN_QC4}
${locUKB3456_QC4}${nl}${locUKB3456_QC4}
${locUKB_ESDPW_QC4}${nl}${locUKB_ESDPW_QC4}
${locUKB_CCPD_QC4}${nl}${locUKB_CCPD_QC4}
${locUKB20161_QC4}${nl}${locUKB20161_QC4}
${locICC_QC4}${nl}${locICC_QC4}
${locUKB_caffeine_QC4}${nl}${locUKB_caffeine_QC4}
_EOF_
) > ${locMR}/filePath_folder_clumped-SNPs_leftJoin_QCed-GWAS-files # 22 lines

# Duplicate each line of QCed GWAS file paths and then combine the two files (22 lines)
## The delimiter is "\t" in the output file
paste -d"\t" ${locMR}/filePath_clumped-SNPs <(awk '{for(i=0;i<2;i++)print}' ${locMR}/filePath_QCed-GWAS-files) <(awk '{for(i=0;i<2;i++)print}' ${locMR}/field-number-SNP_delimiters_QCed-GWAS-files) ${locMR}/filePath_folder_clumped-SNPs_leftJoin_QCed-GWAS-files > ${locMR}/filePath_clumped-SNPs_QCed-GWAS-files

#----------------------------------------------------------------------------------------------------------------#
#------------------------Archive existing files in the output folders--------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
# IFS=$'\n';
# ## Number of iterations: 20
# count=0;
# for line in `cat ${locMR}/filePath_clumped-SNPs_QCed-GWAS-files`; do
# count=$((${count}+1)); 
# echo "========================================== iteration ${count} ============================================"; 
# output_folder=`echo $line | awk '{ print $3}'`;
# echo "output_folder=$output_folder";
# cd $output_folder;
# mkdir -p $output_folder/archive_on_20190404
# mv !(archive_on_20190404) $output_folder/archive_on_20190404/
# done

#----------------------------------------------------------------------------------------------------------------#
#------------------------LeftJoin file1 and file2 ---------------------------------------------------------------#
#--------------------------------- file1 : Clumped SNPs files (expect SNP in field 3)----------------------------#
#--------------------------------- file2 : QCed GSCAN GWAS file of the same trait as file1-----------------------#
#----------------------------------------------------------------------------------------------------------------#

# Loop through every line of a two-column file. In each iteration, leftJoin $fileLeft & $fileRight
## reference URL: https://slaptijack.com/programming/two-column-for-loop-in-bash.html 
## The key here is to set the internal field separator ($IFS) to $'\n' so that the for loop interates on lines rather than words
IFS=$'\n';
## Number of iterations: 22
count=0;
#for line in `sed -n '21,22p' ${locMR}/filePath_clumped-SNPs_QCed-GWAS-files`;do
for line in `cat ${locMR}/filePath_clumped-SNPs_QCed-GWAS-files`; do 
count=$((${count}+1)); 
echo "========================================== iteration ${count} ============================================"; 
filePath_left_table=`echo $line | awk '{ print $1}'`;
filePath_right_table=`echo $line | awk '{ print $2}'`;
SNP_field_numb_fileRight=`echo $line | awk '{ print $3}'`;
delimiter_fileRight=`echo $line | awk '{ print $4}'`;
output_folder=`echo $line | awk '{ print $5}'`;
fileName_fileRight_part1=`basename $filePath_right_table | cut -d'.' -f1`;
# Make output file name
fileName_fileLeft=`basename $filePath_left_table`;
# Remove .clumped using parameter expansion
fileName_fileLeft_extRm=`echo "${fileName_fileLeft%.*}"`;
outputFileName="GWAS_from-clumped-SNPs_${fileName_fileLeft_extRm}";
echo "filePath_left_table=$filePath_left_table"; 
echo "filePath_right_table=$filePath_right_table"; 
echo "SNP_field_numb_fileRight=$SNP_field_numb_fileRight";
echo "output_folder=$output_folder"
echo "fileName_fileRight_part1=$fileName_fileRight_part1"
echo "fileName_fileLeft_extRm=$fileName_fileLeft_extRm"
echo "outputFileName=$outputFileName"
echo "awk '{print \$3}' $filePath_left_table | sed '/^ *\$/d' > ${output_folder}/clumped-SNPs_${fileName_fileRight_part1}"
echo "$locFileJoiner -quiet ${output_folder}/clumped-SNPs_${fileName_fileRight_part1},headerline,key=1,outfields=1 ${filePath_right_table},headerline,sep=${delimiter_fileRight},key=${SNP_field_numb_fileRight},outfields=1- > ${output_folder}/$outputFileName;";
# Extract the SNP column from clumped files
## sed '/^ *$/d' here deletes blank lines at the end of a file
awk '{print $3}' $filePath_left_table | sed '/^ *$/d' > ${output_folder}/clumped-SNPs_${fileName_fileRight_part1} ;
$locFileJoiner -quiet ${output_folder}/clumped-SNPs_${fileName_fileRight_part1},headerline,key=1,outfields=1 ${filePath_right_table},headerline,sep=${delimiter_fileRight},key=${SNP_field_numb_fileRight},outfields=1- > ${output_folder}/$outputFileName;
done

##---------------------------------This is the end of this file-------------------------------------------------##