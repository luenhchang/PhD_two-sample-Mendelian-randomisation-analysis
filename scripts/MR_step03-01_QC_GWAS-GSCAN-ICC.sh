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
## 20190408	exported ${locQC3}/${shortFileName}.ambiguousSNPRemoved. Inconsistent number of header columns (13) and data columns (14) found. Kept first 12 columns
## 20181119	QC GSCAN GWAS files (sample excluded ICC, UKB, 23andme, QIMR BLTS)
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input	${locGSCAN}/ai_noICC.txt
# Input	${locGSCAN}/cpd_noICC.txt
# Input	${locGSCAN}/dpw_noICC.txt
# Input	${locGSCAN}/sc_noICC.txt
# Input	${locGSCAN}/si_noICC.txt
# Input ${locMR}/Cannabis_ICC_UKB.txt
# Outpu ${locQC3}/ai_noICC.ambiguSNPRemoved
# Outpu ${locQC3}/cpd_noICC.ambiguSNPRemoved
# Outpu ${locQC3}/dpw_noICC.ambiguSNPRemoved
# Outpu ${locQC3}/sc_noICC.ambiguSNPRemoved
# Outpu ${locQC3}/si_noICC.ambiguSNPRemoved
# Outpu ${locMR}/Cannabis_ICC_UKB_small.txt
#---------------------------------------------------------------------------------------------------------------
## Locations of main folders
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";
locGSCAN="${homeDir}/LabData/Lab_NickM/lunC/GSCAN";
locICC="${homeDir}/LabData/Lab_NickM/lunC/international_cannabis_consortium";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";
locGSCAN="$locMR/noICC_results";
locRaw=$locGSCAN/QC0_rawdata;
locQC1=$locGSCAN/QC1_find_allOccurencesOfDuplicatedSNPs;
locQC2=$locGSCAN/QC2_remove_duplicatedSNPs;
locQC3=$locGSCAN/QC3_remove_ambiguousSNPs_indel;

mkdir -p $locRaw $locQC1 $locQC2 $locQC3;

# Move raw data to a folder
cp -n $locGSCAN/*_noICC.txt $locRaw

# Store file paths of input files as one file
realpath $locRaw/*_noICC.txt > $locRaw/filePath

#-----------------------------------------------------------------------------------------------------
# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
#-----------------------------------------------------------------------------------------------------
field_RSNum=3; # filed position of SNP Rs ID
for filePath in `cat $locRaw/filePath`; do 
fileName=`basename $filePath`; 
echo $fileName; 
awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $locQC1/${fileName}.allOccurrenceDup; 
done;

#-----------------------------------------------------------------------------------------------------
# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
#-----------------------------------------------------------------------------------------------------
## output files are headerless
field_RSNum=3; # filed position of SNP Rs ID
for filePath in `cat $locRaw/filePath`; do 
fileName=`basename $filePath`; 
echo $fileName;
awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $locQC2/${fileName}.allOccurrenceDupRemoved;
done

## Store file paths as a file
realpath $locQC2/*.allOccurrenceDupRemoved > $locQC2/filePath

#-----------------------------------------------------------------------------------------------------
# Further remove non SNPid (rs number) in the SNP field from files from previous step
#-----------------------------------------------------------------------------------------------------
# Then remove ambiguous SNPs and insertions, deletions
## Number of iterations: 5
for filePath in `cat $locQC2/filePath`;do 
shortFileName=`basename $filePath | cut -d"." -f1`; 
echo "filePath=$filePath"
echo "shortFileName=$shortFileName" 
# Keep the first 12 columns of the data. The header has 13 columns but data has 14 columns.
awk 'BEGIN {FS="\t";OFS="\t"; print "CHROM","POS","RSID","REF","ALT","STAT","PVALUE","BETA","SE","N","EFFECTIVE_N","Number_of_Studies"} {if( $3 ~/rs/) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $filePath | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${locQC3}/${shortFileName}.ambiguousSNPRemoved;
done

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Select wanted columns from cannabis GWAS. In this data set, there are 8% of data with 10 columns and 92% with 11 columns, as checked in the following code
#------------------------------------------------------------------------------------------------------------------------------------------------------
## count line number and column numbers
cd $locMR
awk '{print NF}' Cannabis_ICC_UKB.txt | sort| uniq -c
# 900134 10
#10635459 11

## Missing values occur in $9 and $10. Here subset important column 1-7 and 11 (SNP Allele1 Allele2 MAF Effect StdErr P N) for subsequent MR analysis
cut -d" " -f1,2,3,4,5,6,7 Cannabis_ICC_UKB.txt > temp1
awk '{print $NF}' Cannabis_ICC_UKB.txt > temp2 # $NF is $10 or $11, depending on whether awk reads the line with missing column. Don't use $11 in the awk, it is different from $NF, which can change line by line
paste -d" " temp1 temp2 > temp3
grep -E "SNP|rs" temp3 > temp4

# Check if lines are still have inconsistent column numbers
awk '{print NF}' temp4 | sort | uniq -c

wc -l temp1 temp2 temp3 temp4
  # 11535593 temp1
  # 11535593 temp2
  # 11535593 temp3
  # 11534877 temp4

# rename file
cat temp4 > Cannabis_ICC_UKB_small.txt

cp -n ${locScripts}/MR_step03_jobSubmit_QC-GWAS.sh ${locScripts}/MR_step04_LD-based-SNP-clumping.sh
cp -n ${locScripts}/MR_step03-01_QC_GWAS-GSCAN-ICC.sh ${locScripts}/MR_step03-02_QC_GWAS-UKB.sh
cp -n ${locScripts}/MR_step03-01_QC_GWAS-GSCAN-ICC.sh ${locScripts}/MR_step03-04_compare-sample-sizes_QC-GWAS-GSCAN_different-samples.sh

##---------------------------------This is the end of this file-------------------------------------------------##