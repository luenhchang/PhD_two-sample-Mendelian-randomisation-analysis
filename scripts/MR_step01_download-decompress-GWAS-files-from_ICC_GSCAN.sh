#!/bin/bash
## file name: MR_step01_download-decompress-GWAS-files-from_ICC_GSCAN.sh
## modified from: ${homeDir}/scripts/PRS_UKB_201711/PRS_UKB_201711_step01_download-GSCAN-GWAS_update-GWAS-info.sh
## date created: 20180628
## purpose: (1) unzip GWAS summary statistics file manually downloaded from international cannabis consortium via SURFDrive
##	(2) unzip GWAS summary statistics files from GSCAN via scp
## Run dependency: 

## How to run this file: 
## Time 	Change
##---------------------------------------------------------------------------------------------------------------------------------------------------
## 20180627 Downloaded cannabis GWAS file from https://surfdrive.surf.nl/files/index.php/s/T282wqYPyJcBQ8Y/download to laptop with password: 532TZma+ (valid until 20180708) and laptop to remote directory ${locICC} 
##---------------------------------------------------------------------------------------------------------------------------------------------------

# Type	Files																Size(MB)
#------------------------------------------------------------------------------------------------------------------------
# Input	${locICC}/Cannabis_ICC_UKB.txt.gz									217.7
# Input	${locGSCAN}/data/shared_results/no23andMe/no23andMe_results.tar.gz	1977
# Input ${locGSCAN}/data/shared_results/noICC/noICC_results.tar.gz			2033
# Outpu	$locMR/no23andMe_results.tar.gz										1977
# Outpu	$locMR/no23andMe_results/ai_no23andme.txt.gz						400
# Outpu	$locMR/no23andMe_results/ai_no23andme.txt							1655
# Outpu	$locMR/no23andMe_results/cpd_no23andme.txt.gz						401
# Outpu	$locMR/no23andMe_results/cpd_no23andme.txt							1658
# Outpu	$locMR/no23andMe_results/dpw_no23andme.txt.gz						402
# Outpu	$locMR/no23andMe_results/dpw_no23andme.txt							1645
# Outpu	$locMR/no23andMe_results/sc_no23andme.txt.gz						375
# Outpu	$locMR/no23andMe_results/sc_no23andme.txt							1568
# Outpu	$locMR/no23andMe_results/si_no23andme.txt.gz						392
# Outpu	$locMR/no23andMe_results/si_no23andme.txt							1623
# Outpu	$locMR/Cannabis_ICC_UKB.txt.gz 										217.7 
# Outpu	$locMR/Cannabis_ICC_UKB.txt 										695
# Outpu	$locMR/noICC_results/ai_noICC.txt									1974
# Outpu	$locMR/noICC_results/cpd_noICC.txt									1966
# Outpu	$locMR/noICC_results/dpw_noICC.txt									2076
# Outpu	$locMR/noICC_results/sc_noICC.txt									1940
# Outpu	$locMR/noICC_results/si_noICC.txt									2007
# Outpu	$locMR/noICC_results/README.txt
# Outpu	$locMR/noICC_results/sc_noICC.txt.gz.tbi							
# Outpu	$locMR/noICC_results/si_noICC.txt.gz.tbi
# Outpu	$locMR/noICC_results/ai_noICC.txt.gz.tbi
# Outpu	$locMR/noICC_results/cpd_noICC.txt.gz.tbi
# Outpu	$locMR/noICC_results/dpw_noICC.txt.gz.tbi
#------------------------------------------------------------------------------------------------------------------

## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";
locGSCAN="${homeDir}/LabData/Lab_NickM/lunC/GSCAN";
locICC="${homeDir}/LabData/Lab_NickM/lunC/international_cannabis_consortium";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";
locMR_MS_SU="${workingDir}/MR_MS_SU_201901/data";
mkdir -p ${locMR} ${locMR_MS_SU};

#---------------------------------------------------------------------------------------------------#
#-------------------------------Download GWAS files from GSCAN server to QIMR LabData folder--------#
#---------------------------------------------------------------------------------------------------#
## Note /mnt/backedup/home/lunC/LabData-fixed will be replaced by /mnt/backedup/home/lunC/LabData
## Username: sftp-gscan
## Host: share.sph.umich.edu
## password: T3g3j8cA
## filePath : ${locGSCAN}/data/shared_results/no23andMe/no23andMe_results.tar.gz
## file size: 1977 MB
cp -n ${locGSCAN}/data/shared_results/no23andMe/no23andMe_results.tar.gz $locMR/no23andMe_results.tar.gz;

# Download GSCAN meta-analysis files on a sample independent of ICC and UKB
## Create a same named folder as the destination of the download 
mkdir -p ${locGSCAN}/data/shared_results/noICC ;

## Map network drive LabData by typing LabData

## Username: sftp-gscan
## Host: share.sph.umich.edu
## password: Q2XK0wohiutCmz7n
## File directory: /data/shared_results/noICC/noICC_results.tar.gz
scp sftp-gscan@share.sph.umich.edu:/data/shared_results/noICC/noICC_results.tar.gz ${locGSCAN}/data/shared_results/noICC/; 
# noICC_results.tar.gz                                                                                                 100% 1939MB   1.8MB/s   17:39

cp -n ${locGSCAN}/data/shared_results/noICC/noICC_results.tar.gz ${locMR}/noICC_results.tar.gz

#---------------------------------------------------------------------------------------------------#
#-------------------------------Give permission to a specific user----------------------------------#
#---------------------------------------------------------------------------------------------------#

## Client username: yuanZ
## File to give permision: /mnt/backedup/home/lunC/LabData/Lab_NickM/lunC/GSCAN/data/shared_results/no23andMe/no23andMe_results.tar.gz
cp -n ${locGSCAN}/data/shared_results/no23andMe/no23andMe_results.tar.gz /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/no23andMe_results.tar.gz

# Allow users to read the file only
chmod 644 /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_GSCAN/no23andMe_results.tar.gz

#---------------------------------------------------------------------------------------------------#
#------------Manually downloaded cannabis GWAS file from ICC SURFDrive to a local folder------------#
#------------Copy tar.gz file to Chang's working folder---------------------------------------------#
#---------------------------------------------------------------------------------------------------#
## file path local: D:\Now\library_genetics_epidemiology_GWAS_largeFiles\international_cannabis_consortium\Cannabis_ICC_UKB.txt.gz
## file path remote: $locICC/Cannabis_ICC_UKB.txt.gz
## file path under working folder: $locMR/Cannabis_ICC_UKB.txt.gz
cp -n $locICC/Cannabis_ICC_UKB.txt.gz $locMR/Cannabis_ICC_UKB.txt.gz;

#---------------------------------------------------------------------------------------------------------------#
# --------------------------------------Decompress tar.gz and txt.gz files
#---------------------------------------------------------------------------------------------------------------#
# Decompress tar.gz files
## file.tar.gz will become file.tar. .tar is a directory. You will tar the directory to see what is inside
dir_source="$locICC";
dir_destin="$locMR";
cd ${dir_destin}/;

#gunzip no23andMe_results.tar.gz; # no23andMe_results.tar
gunzip noICC_results.tar.gz; # noICC_results.tar

# Untar the tar archive file
## Go to the file directory so the untar file will be generated there. Otherwise, it will be in your home folder
#tar -xvf no23andMe_results.tar;
# no23andMe_results/
# no23andMe_results/ai_no23andme.txt.gz
# no23andMe_results/cpd_no23andme.txt.gz
# no23andMe_results/dpw_no23andme.txt.gz
# no23andMe_results/sc_no23andme.txt.gz
# no23andMe_results/si_no23andme.txt.gz
# no23andMe_results/sc_no23andme.txt.gz.tbi
# no23andMe_results/si_no23andme.txt.gz.tbi
# no23andMe_results/dpw_no23andme.txt.gz.tbi
# no23andMe_results/cpd_no23andme.txt.gz.tbi
# no23andMe_results/ai_no23andme.txt.gz.tbi
# no23andMe_results/README.txt
tar -xvf noICC_results.tar
# noICC_results/
# noICC_results/ai_noICC.txt.gz
# noICC_results/cpd_noICC.txt.gz
# noICC_results/dpw_noICC.txt.gz
# noICC_results/sc_noICC.txt.gz
# noICC_results/si_noICC.txt.gz
# noICC_results/ai_noICC.txt.gz.tbi
# noICC_results/cpd_noICC.txt.gz.tbi
# noICC_results/sc_noICC.txt.gz.tbi
# noICC_results/dpw_noICC.txt.gz.tbi
# noICC_results/si_noICC.txt.gz.tbi
# noICC_results/README.txt

# Decompress txt.gz files
## txt.gz files become txt files
gunzip Cannabis_ICC_UKB.txt.gz;  

## Decompress *_noICC.txt.gz files within the folder
#cd ${dir_destin}/no23andMe_results;
cd ${dir_destin}/noICC_results

#for file in `ls *_no23andme.txt.gz`;do echo $file; gunzip $file; done;
for file in $(ls *_noICC.txt.gz); do 
echo $file;
gunzip $file;
done

# Copy this script for similar jobs
#cp -n ${homeDir}/scripts/MR_step01_download-decompress-GWAS-files-from_ICC_GSCAN.sh
#------------------------------------------------------------------------------------------------------#
#----------------------------------This is the end of this file----------------------------------------#
#------------------------------------------------------------------------------------------------------#