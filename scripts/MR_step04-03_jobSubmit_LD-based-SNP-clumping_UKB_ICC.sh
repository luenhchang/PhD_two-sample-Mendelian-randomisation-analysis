#!/bin/bash

## file name: MR_step04-03_jobSubmit_LD-based-SNP-clumping_UKB_ICC.sh
## old file name: zMR_step04-03_jobSubmit_LD-based-SNP-clumping_UKB.sh
## modified from: MR_step04-02_jobSubmit_LD-based-SNP-clumping_GSCAN.sh
## date created: 20190403
## purpose: select independent SNPs
## Run dependency: 
## How to run this script: . ${locScripts}/MR_step04-03_jobSubmit_LD-based-SNP-clumping_UKB.sh > ${locHistory}/jobSubmitted_20190403_LD-based-clumping_UKB-GWASs

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20190812	qsub jobs 8593545-8593546 (2 jobs, UKB caffeine only)
## 20190404	qsub jobs 7277912-7277916 (2 jobs, ICC only)
## 20190403	qsub jobs 7275949-7275956 (8 jobs, UKB only)
## 20180823	qsub jobs (4*3)
## 20180822	qsub jobs 6092089-6092092 (4 jobs) 6092126-6092129 (4 jobs)
## 20180804 qsub jobs 6048942-6048945 (4 jobs)
## --clump-p1 5e-8 or 1e-6 [index variant p-value threshold] This is the p value threshold for choosing the best SNP per haplotype
## --clump-p2 1e-6 [clumped variant p-value threshold] Once the best SNP is chosen, other SNPs within the same haplotype with p value < p2 and r2 0.01 are chosen (i.e. clumped)
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input	${locInput}/LD-SNP-clumping-criteria.txt

# Input $locUKB3456_QC3/QCed-GWAS-UKB3456_headed
# Input $locUKB_ESDPW_QC3/QCed-GWAS-UKB-ESDPW_headed
# Input $locUKB_CCPD_QC3/QCed-GWAS-UKB-CCPD_headed
# Input $locUKB20161_QC3/QCed-GWAS-UKB-PYOS_headed
# Input $locICC/Cannabis_ICC_UKB_small.txt
# Input ${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed

# Outpu	$locUKB3456_input/SNP_PValue_QCed-GWAS-UKB3456_headed
# Outpu	$locUKB_ESDPW_input/SNP_PValue_QCed-GWAS-UKB-ESDPW_headed
# Outpu	$locUKB_CCPD_input/SNP_PValue_QCed-GWAS-UKB-CCPD_headed
# Outpu	$locUKB20161_input/SNP_PValue_QCed-GWAS-UKB-PYOS_headed
# Outpu	$locICC_input/SNP_PValue_Cannabis_ICC_UKB_small
# Outpu	${locUKB_caffeine_input}/SNP_PValue_QCed-GWAS-UKB-caffeine-consumed-per-day_headed

# Outpu	$locMR/filePath_QCed-UKB-ICC-GWAS_SNP-PValue-subset.txt
# Outpu ${locUKB3456_output}/GWAS-UKB3456_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Outpu $locUKB_ESDPW_output/GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Outpu $locUKB_CCPD_output/GWAS-UKB-CCPD_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Outpu $locUKB20161_output/GWAS-UKB-PYOS_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Outpu $locICC_output/GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
# Outpu $locUKB_caffeine_output/GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-*_p2-*.clumped (2 files)
#---------------------------------------------------------------------------------------------------------------

#-----------------------------------------------
# Folder locations under my home 
#-----------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
jobScriptFilePath="${locScripts}/MR_step04-01-01_jobScript_LD-based-SNP-clumping.sh";
locHistory="${homeDir}/history";

#-----------------------------------------------
# Folder locations under my working
#-----------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";

# Folder locations related to ICC cannabis initiation
locICC="$locMR/ICC-cannabis-ever";
locICC_LDBasedClumping="$locICC/LD-based_SNP_clumping";
locICC_input=$locICC_LDBasedClumping/input;
locICC_output=$locICC_LDBasedClumping/output;
locICC_PBS_output=$locICC_LDBasedClumping/pbs_output;
mkdir -p $locICC_LDBasedClumping $locICC_input $locICC_output $locICC_PBS_output ;

# Folder locations related to GSCAN
locGSCAN="$locMR/noICC_results";
locLDBasedClumping=$locGSCAN/LD-based_SNP_clumping;
locInput=$locLDBasedClumping/input;

# Folder locations related to UKB3456 exluding people with ever using cannabis data
locUKB3456=$locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis
locUKB3456_QC3=$locUKB3456/QC3_remove_ambiguousSNPs_indel;

locUKB3456_LDBasedClumping=$locUKB3456/LD-based_SNP_clumping;
locUKB3456_input=$locUKB3456_LDBasedClumping/input;
locUKB3456_output=$locUKB3456_LDBasedClumping/output;
locUKB3456_PBS_output=$locUKB3456_LDBasedClumping/pbs_output;

# Folder locations related to UKB ESDPW exluding people with ever using cannabis data
locUKB_ESDPW=$locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis;
locUKB_ESDPW_QC3=$locUKB_ESDPW/QC3_remove_ambiguousSNPs_indel;
locUKB_ESDPW_LDBasedClumping=$locUKB_ESDPW/LD-based_SNP_clumping;
locUKB_ESDPW_input=$locUKB_ESDPW_LDBasedClumping/input;
locUKB_ESDPW_output=$locUKB_ESDPW_LDBasedClumping/output;
locUKB_ESDPW_PBS_output=$locUKB_ESDPW_LDBasedClumping/pbs_output;

# Folder locations related to UKB CCPD exluding people with ever using cannabis data
locUKB_CCPD=$locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis;
locUKB_CCPD_QC3=$locUKB_CCPD/QC3_remove_ambiguousSNPs_indel;
locUKB_CCPD_LDBasedClumping=$locUKB_CCPD/LD-based_SNP_clumping;
locUKB_CCPD_input=$locUKB_CCPD_LDBasedClumping/input;
locUKB_CCPD_output=$locUKB_CCPD_LDBasedClumping/output;
locUKB_CCPD_PBS_output=$locUKB_CCPD_LDBasedClumping/pbs_output;

# Folder locations related to UKB 20161 PYOS exluding people with ever using cannabis data
locUKB20161=$locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis;
locUKB20161_QC3=$locUKB20161/QC3_remove_ambiguousSNPs_indel;
locUKB20161_LDBasedClumping=$locUKB20161/LD-based_SNP_clumping;
locUKB20161_input=$locUKB20161_LDBasedClumping/input;
locUKB20161_output=$locUKB20161_LDBasedClumping/output;
locUKB20161_PBS_output=$locUKB20161_LDBasedClumping/pbs_output;

# Folder locations related to UKB estimated caffeine consumed per day
locUKB_caffeine=${locMR}/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis;
locUKB_caffeine_raw=${locUKB_caffeine}/QC0_rawdata;
locUKB_caffeine_QC1=${locUKB_caffeine}/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB_caffeine_QC2=${locUKB_caffeine}/QC2_remove_duplicatedSNPs;
locUKB_caffeine_QC3=${locUKB_caffeine}/QC3_remove_ambiguousSNPs_indel;
locUKB_caffeine_LDBasedClumping=${locUKB_caffeine}/LD-based_SNP_clumping ;
locUKB_caffeine_input=${locUKB_caffeine_LDBasedClumping}/input;
locUKB_caffeine_output=${locUKB_caffeine_LDBasedClumping}/output;
locUKB_caffeine_PBS_output=${locUKB_caffeine_LDBasedClumping}/pbs_output;

mkdir -p $locUKB3456_LDBasedClumping $locUKB3456_input $locUKB3456_output $locUKB3456_PBS_output;
mkdir -p $locUKB_ESDPW_QC3 $locUKB_ESDPW_LDBasedClumping $locUKB_ESDPW_input $locUKB_ESDPW_output $locUKB_ESDPW_PBS_output;   
mkdir -p $locUKB_CCPD_QC3 $locUKB_CCPD_LDBasedClumping $locUKB_CCPD_input $locUKB_CCPD_output $locUKB_CCPD_PBS_output;
mkdir -p $locUKB20161_QC3 $locUKB20161_LDBasedClumping $locUKB20161_input $locUKB20161_output $locUKB20161_PBS_output;
mkdir -p ${locUKB_caffeine_LDBasedClumping} ${locUKB_caffeine_input} ${locUKB_caffeine_output} ${locUKB_caffeine_PBS_output}; 

bfile="/mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/LD_reference/LD_ref_201803_rsOnly";
email_address="luenhchang@gmail.com";

#------------------------------------------------------------------------------------------------------------------------------#
#---Subset SNP and P columns for LD clumping from input files
## Only do this once
#------------------------------------------------------------------------------------------------------------------------------#
# awk 'BEGIN {print "SNP P"} (NR !=1) {print $1, $16}' $locUKB3456_QC3/QCed-GWAS-UKB3456_headed > $locUKB3456_input/SNP_PValue_QCed-GWAS-UKB3456_headed;
# awk 'BEGIN {print "SNP P"} (NR !=1) {print $1, $16}' $locUKB_ESDPW_QC3/QCed-GWAS-UKB-ESDPW_headed > $locUKB_ESDPW_input/SNP_PValue_QCed-GWAS-UKB-ESDPW_headed;
# awk 'BEGIN {print "SNP P"} (NR !=1) {print $1, $16}' $locUKB_CCPD_QC3/QCed-GWAS-UKB-CCPD_headed > $locUKB_CCPD_input/SNP_PValue_QCed-GWAS-UKB-CCPD_headed;
# awk 'BEGIN {print "SNP P"} (NR !=1) {print $1, $16}' $locUKB20161_QC3/QCed-GWAS-UKB-PYOS_headed > $locUKB20161_input/SNP_PValue_QCed-GWAS-UKB-PYOS_headed;
# awk 'BEGIN {print "SNP P"} (NR !=1) {print $1, $7}' $locICC/Cannabis_ICC_UKB_small.txt > $locICC_input/SNP_PValue_Cannabis_ICC_UKB_small
# awk 'BEGIN {print "SNP P"} (NR !=1) {print $1, $16}' ${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed > ${locUKB_caffeine_input}/SNP_PValue_QCed-GWAS-UKB-caffeine-consumed-per-day_headed;

# Write 3 columns into a file for iterations 
## Column 1: file paths of QCed GWAS files
## Column 2: Subset files (SNP, P values) of the QCed GWAS files
## Column 3: output folder path for clumped files
## Column 4: User-named prefix for output file
(cat <<- _EOF_
GWAS_file_path SNP_PValue_file_path outputFolderPath outputFileNamePrefix
$locUKB3456_QC3/QCed-GWAS-UKB3456_headed $locUKB3456_input/SNP_PValue_QCed-GWAS-UKB3456_headed $locUKB3456_output GWAS-UKB3456
$locUKB_ESDPW_QC3/QCed-GWAS-UKB-ESDPW_headed $locUKB_ESDPW_input/SNP_PValue_QCed-GWAS-UKB-ESDPW_headed $locUKB_ESDPW_output GWAS-UKB-ESDPW
$locUKB_CCPD_QC3/QCed-GWAS-UKB-CCPD_headed $locUKB_CCPD_input/SNP_PValue_QCed-GWAS-UKB-CCPD_headed $locUKB_CCPD_output GWAS-UKB-CCPD
$locUKB20161_QC3/QCed-GWAS-UKB-PYOS_headed $locUKB20161_input/SNP_PValue_QCed-GWAS-UKB-PYOS_headed $locUKB20161_output GWAS-UKB-PYOS
$locICC/Cannabis_ICC_UKB_small.txt $locICC_input/SNP_PValue_Cannabis_ICC_UKB_small $locICC_output GWAS-ICC-CI
${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed ${locUKB_caffeine_input}/SNP_PValue_QCed-GWAS-UKB-caffeine-consumed-per-day_headed ${locUKB_caffeine_output} GWAS-UKB-caffeine
_EOF_
) > $locMR/filePath_QCed-UKB-ICC-GWAS_SNP-PValue-subset.txt

#------------------------------------------------------------------------------------------------------------------------------#
# Archive existing files to the archive folder
## ONLY do this once
#------------------------------------------------------------------------------------------------------------------------------#
# IFS=$'\n';
# for line in `tail -n+2 $locMR/filePath_QCed-UKB-GWAS_SNP-PValue-subset.txt`; do
# outputFolderPath=$(echo $line | cut -d" " -f3);
# pbs_output_folderPath=$(dirname $outputFolderPath)/pbs_output;
# echo "outputFolderPath=${outputFolderPath}";
# echo "pbs_output_folderPath=${pbs_output_folderPath}";
# mkdir -p ${outputFolderPath}/archive ${pbs_output_folderPath}/archive;
# # Archive all files in ${outputFolderPath} except the archive folder
# cd ${outputFolderPath};
# mv !(archive) ${outputFolderPath}/archive/ ;
# # Archive all files in $pbs_output_folderPath except the archive folder
# cd ${pbs_output_folderPath};
# mv !(archive) ${pbs_output_folderPath}/archive/ ;
# done

#------------------------------------------------------------------------------------------------------------------------------#
# Perform LD clumping for GWASs from UKB, ICC
## Loop thru each line of 2 space-separated files, skipping the first row of header text
## The input files are SNP and P value files
#------------------------------------------------------------------------------------------------------------------------------#
## Number of iteration: 6 (files) * 2 (sets of clumping criteria)
IFS=$'\n';
count=0;
for line1 in `tail -n+2 $locMR/filePath_QCed-UKB-ICC-GWAS_SNP-PValue-subset.txt`; do
#for line1 in `tail -n 1 $locMR/filePath_QCed-UKB-ICC-GWAS_SNP-PValue-subset.txt`; do
GWAS_file_path=$(echo $line1 | cut -d" " -f1);
SNP_PValue_file_path=$(echo $line1 | cut -d" " -f2);
outputFolderPath=$(echo $line1 | cut -d" " -f3);
pbs_output_folderPath=$(dirname $outputFolderPath)/pbs_output
outputFileNamePrefix=$(echo $line1 | cut -d" " -f4);
for line2 in `tail -n+2 ${locInput}/LD-SNP-clumping-criteria.txt`; do
count=$((${count}+1)); 
jobName="job${count}"
clumping_r2_threshold=$(echo $line2 | cut -d" " -f1);
LDWindow=$(echo $line2 | cut -d" " -f2);
p_threshold_1=$(echo $line2 | cut -d" " -f3);
p_threshold_2=$(echo $line2 | cut -d" " -f4);
outputFileName="${outputFileNamePrefix}_LDWindow-kb-${LDWindow}_R2-${clumping_r2_threshold}_p1-${p_threshold_1}_p2-${p_threshold_2}";
echo "========================================== iteration ${count} ============================================"; 
echo "GWAS_file_path=$GWAS_file_path";
echo "SNP_PValue_file_path=$SNP_PValue_file_path";
echo "outputFolderPath=$outputFolderPath";
echo "pbs_output_folderPath=$pbs_output_folderPath";
echo "outputFileNamePrefix=$outputFileNamePrefix";
echo "clumping_r2_threshold=$clumping_r2_threshold";
echo "LDWindow=$LDWindow";
echo "p_threshold_1=$p_threshold_1";
echo "p_threshold_2=$p_threshold_2";
echo "outputFileName=${outputFileName}";
echo "qsub -N ${jobName} -m bea -M $email_address -v v_bfilePath=${bfile},v_clumpFilePath=${SNP_PValue_file_path},v_clumping_r2_threshold=${clumping_r2_threshold},v_clumping_LD_window=${LDWindow},v_p_threshold_1=${p_threshold_1},v_p_threshold_2=${p_threshold_2},v_outputDir=${outputFolderPath},v_outputResultFileName=${outputFileName} -e ${pbs_output_folderPath}/${outputFileName}.pbs.err -o ${pbs_output_folderPath}/${outputFileName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=2500mb ${jobScriptFilePath};"
qsub -N ${jobName} -m bea -M $email_address -v v_bfilePath=${bfile},v_clumpFilePath=${SNP_PValue_file_path},v_clumping_r2_threshold=${clumping_r2_threshold},v_clumping_LD_window=${LDWindow},v_p_threshold_1=${p_threshold_1},v_p_threshold_2=${p_threshold_2},v_outputDir=${outputFolderPath},v_outputResultFileName=${outputFileName} -e ${pbs_output_folderPath}/${outputFileName}.pbs.err -o ${pbs_output_folderPath}/${outputFileName}.pbs.out -l ncpus=1,walltime=01:00:00,mem=2500mb ${jobScriptFilePath};
done
done


##---------------------------------This is the end of this file-------------------------------------------------##
