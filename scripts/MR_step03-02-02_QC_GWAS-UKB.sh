#!/bin/bash
## File path: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step03-02-02_QC_GWAS-UKB.sh 
## old file name: 
## Modified from: MR_step03-01_QC_GWAS-GSCAN-ICC.sh
## date created: 20180804
## purpose: Quality Control GWAS files by removing duplicated SNPs, ambiguous SNPs
## Run dependency: 
## Functions external: copy_files_get_filePaths, remove_multiple_occurrences_SNPs, remove_nonSNPid_ambiguous_SNPs_insertions_deletions
## How to run this script: 

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20190812	Exported ${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed
## 20180823 Generated ${locUKB_CCPD_QC3}/QCed-GWAS-UKB-CCPD_headed, ${locUKB20161_QC3}/QCed-GWAS-UKB-PYOS_headed, $locUKB_ESDPW_QC3/QCed-GWAS-UKB-ESDPW_headed  
## 20180821	Generated $locUKB_ESDPW_QC3/QCed-GWAS-UKB-ESDPW_headed
## 20180820 Generated ${locUKB_CCPD_QC3}/QCed-GWAS-UKB-CCPD_headed
## 20180804	Generated ${locUKB3456_QC3}/QCed-GWAS-UKB3456_headed
##--------------------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input ${ref_UKB_3456_NA20453_GWAS}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_X3456_mean.bgen.assoc
# Input ${ref_UKB_ESDPW_NA20453_GWAS}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_complete_alcohol_unitsweekly.bgen.assoc
# Input ${ref_UKB_20161_NA20453_GWAS}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_merged_pack_years_20161.bgen.assoc
# Input ${ref_UKB_CCPD_NA20453_GWAS}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_all_coffee_cpd.bgen.assoc
# Input ${ref_UKB_caffeine_NA20453_GWAS}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_caffeine.per.day.bgen.assoc

# Outpu ${locUKB3456_QC3}/QCed-GWAS-UKB3456_headed
# Outpu ${locUKB_ESDPW_QC3}/QCed-GWAS-UKB-ESDPW_headed
# Outpu ${locUKB_CCPD_QC3}/QCed-GWAS-UKB-CCPD_headed
# Outpu ${locUKB20161_QC3}/QCed-GWAS-UKB-PYOS_headed
# Outpu ${locUKB_caffeine_QC3}/QCed-GWAS-UKB-caffeine-consumed-per-day_headed
# Outpu ${locUKB20453_QC3}/QCed-GWAS-UKB-ever-taken-cannabis_headed
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------
# Folder locations under my home
#---------------------------------------------------
homeDir="/mnt/backedup/home/lunC";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806";
locHistory="${homeDir}/history";

#---------------------------------------------------
# Folder locations under Stuart's lab
#---------------------------------------------------
referenceUKBGWAS="/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas";

ref_UKB_3456_NA20453_GWAS=${referenceUKBGWAS}/BOLT_LMM/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/BOLT-LMM-ukb3456_IID_NA_in_20453-pheno-X3456_mean;

#locUKB_CCPD_GWAS="${referenceUKBGWAS}/BOLT_LMM/coffee/all_coffee_cpd"; # note this GWAS not used as it didn't exclude cannabis users
ref_UKB_ESDPW_NA20453_GWAS="$referenceUKBGWAS/BOLT_LMM/UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/BOLT-LMM-phenotype-pheno-complete_alcohol_unitsweekly";
ref_UKB_20161_NA20453_GWAS=${referenceUKBGWAS}/BOLT_LMM/UKB20161-pack-years-of-smoking_IID-NA-in-UKB20453-everUsedCannabis/BOLT-LMM-phenotype-pheno-merged_pack_years_20161

ref_UKB_CCPD_NA20453_GWAS="${referenceUKBGWAS}/BOLT_LMM/UKB-cups-of-coffee-per-day_IID-NA-in-UKB20453-everUsedCannabis/BOLT-LMM-phenotype-pheno-all_coffee_cpd";

ref_UKB_caffeine_NA20453_GWAS="${referenceUKBGWAS}/BOLT_LMM/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/BOLT-LMM-phenotype-pheno-caffeine.per.day";

ref_UKB_20453_GWAS="${referenceUKBGWAS}/plink2/UKB20453-ever-taken-cannabis/plink2-ukb20453.phenoUtility.recoded-pheno-X20453_0_0_recoded" ;

#---------------------------------------------------
# Folder locations under my working
#---------------------------------------------------
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";

locUKB3456=$locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis
locUKB3456_raw=$locUKB3456/QC0_rawdata;
locUKB3456_QC1=$locUKB3456/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB3456_QC2=$locUKB3456/QC2_remove_duplicatedSNPs;
locUKB3456_QC3=$locUKB3456/QC3_remove_ambiguousSNPs_indel;

locUKB_CCPD=$locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis;
locUKB_CCPD_raw=$locUKB_CCPD/QC0_rawdata;
locUKB_CCPD_QC1=$locUKB_CCPD/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB_CCPD_QC2=$locUKB_CCPD/QC2_remove_duplicatedSNPs;
locUKB_CCPD_QC3=$locUKB_CCPD/QC3_remove_ambiguousSNPs_indel;

locUKB_caffeine=$locMR/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis;
locUKB_caffeine_raw=$locUKB_caffeine/QC0_rawdata;
locUKB_caffeine_QC1=$locUKB_caffeine/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB_caffeine_QC2=$locUKB_caffeine/QC2_remove_duplicatedSNPs;
locUKB_caffeine_QC3=$locUKB_caffeine/QC3_remove_ambiguousSNPs_indel;

locUKB_ESDPW=$locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis;
locUKB_ESDPW_raw=$locUKB_ESDPW/QC0_rawdata;
locUKB_ESDPW_QC1=$locUKB_ESDPW/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB_ESDPW_QC2=$locUKB_ESDPW/QC2_remove_duplicatedSNPs;
locUKB_ESDPW_QC3=$locUKB_ESDPW/QC3_remove_ambiguousSNPs_indel;

locUKB20161=$locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis
locUKB20161_raw=$locUKB20161/QC0_rawdata;
locUKB20161_QC1=$locUKB20161/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB20161_QC2=$locUKB20161/QC2_remove_duplicatedSNPs;
locUKB20161_QC3=$locUKB20161/QC3_remove_ambiguousSNPs_indel;

locUKB20453=$locMR/UKB20453-ever-taken-cannabis
locUKB20453_raw=$locUKB20453/QC0_rawdata;
locUKB20453_QC1=$locUKB20453/QC1_find_allOccurencesOfDuplicatedSNPs;
locUKB20453_QC2=$locUKB20453/QC2_remove_duplicatedSNPs;
locUKB20453_QC3=$locUKB20453/QC3_remove_ambiguousSNPs_indel;


mkdir -p $locUKB3456 $locUKB3456_raw $locUKB3456_QC1 $locUKB3456_QC2 $locUKB3456_QC3 ${locUKB_CCPD} ${locUKB_CCPD_raw} ${locUKB_CCPD_QC1} ${locUKB_CCPD_QC2} ${locUKB_CCPD_QC3} ${locUKB_ESDPW_raw} ${locUKB_ESDPW_QC1} ${locUKB_ESDPW_QC2} ${locUKB_ESDPW_QC3} ${locUKB_caffeine_raw} ${locUKB20161_raw}  ${locUKB20161_QC1} ${locUKB20161_QC2} ${locUKB20161_QC3} ${locUKB_caffeine_QC1} ${locUKB_caffeine_QC2} ${locUKB_caffeine_QC3};

mkdir -p ${locUKB20453} ${locUKB20453_raw} ${locUKB20453_QC1} ${locUKB20453_QC2} ${locUKB20453_QC3}   

# Move raw GWAS data to your working folder
cp -n $ref_UKB_3456_NA20453_GWAS/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_X3456_mean.bgen.assoc $locUKB3456_raw;
cp -n $ref_UKB_ESDPW_NA20453_GWAS/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_complete_alcohol_unitsweekly.bgen.assoc $locUKB_ESDPW_raw;
cp -n $ref_UKB_20161_NA20453_GWAS/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_merged_pack_years_20161.bgen.assoc $locUKB20161_raw;
cp -n $ref_UKB_CCPD_NA20453_GWAS/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_all_coffee_cpd.bgen.assoc $locUKB_CCPD_raw;
cp -n ${ref_UKB_20453_GWAS}/X20453_0_0_recoded_plink-ukb_imp_chr{1..22}_v3.plink2.output.X20453_0_0_recoded_plink.glm.logistic ${locUKB20453_raw};

# Store file paths of input GWAS files as one file
realpath $locUKB3456_raw/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_X3456_mean.bgen.assoc > $locUKB3456_raw/filePath;
realpath $locUKB_ESDPW_raw/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_complete_alcohol_unitsweekly.bgen.assoc > $locUKB_ESDPW_raw/filePath;
realpath $locUKB20161_raw/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_merged_pack_years_20161.bgen.assoc > $locUKB20161_raw/filePath;
realpath $locUKB_CCPD_raw/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_all_coffee_cpd.bgen.assoc > $locUKB_CCPD_raw/filePath;
realpath ${locUKB20453_raw}/X20453_0_0_recoded_plink-ukb_imp_chr{1..22}_v3.plink2.output.X20453_0_0_recoded_plink.glm.logistic > ${locUKB20453_raw}/filePath;

#---------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------- QC UKB GWAS 3456 
#----------------------------------------------------------exclude people with data in ever using cannabis 
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB3456_raw/filePath`; do 
	fileName=`basename $filePath`; 
	echo "fileName=$fileName"; 
	awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $locUKB3456_QC1/${fileName}.allOccurrenceDup; 
done;

# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
## output files are headerless
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB3456_raw/filePath`; do 
	fileName=`basename $filePath`; 
	echo "fileName=$fileName";
	awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $locUKB3456_QC2/${fileName}.allOccurrenceDupRemoved;
done

## Store file paths as a file
realpath $locUKB3456_QC2/*.allOccurrenceDupRemoved > $locUKB3456_QC2/filePath

# Further remove non SNPid (rs number) in the SNP field from files from previous step
# Then remove ambiguous SNPs and insertions, deletions
for filePath in `cat $locUKB3456_QC2/filePath`;do 
	fileName=`basename $filePath`; 
	echo "fileName=$fileName";
	awk 'BEGIN {FS="\t";OFS="\t"; print "SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","P_BOLT_LMM" } {if( $1 ~/rs/) print $0}' $filePath | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${locUKB3456_QC3}/${fileName}.ambiguousSNPRemoved;
done

# Combine 22 autosome GWAS files as a single file
## Save file paths in a file
realpath ${locUKB3456_QC3}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_X3456_mean.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved > ${locUKB3456_QC3}/filePath

## Vertically combine all 22 files excluding their header rows
cat ${locUKB3456_QC3}/filePath | xargs awk '(FNR>1){print $0}'> ${locUKB3456_QC3}/QCed-GWAS-UKB3456

##Write header back to the file using in-place computing
cd $locUKB3456_QC3
file1=revised_bolt_imputed_ukb_imp_chr1_v3_X3456_mean.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved
file2=QCed-GWAS-UKB3456
cat <(head -1 $file1 | sed 's/P_BOLT_LMM/PVALUE/g') <(cat $file2) > QCed-GWAS-UKB3456_headed

#---------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ QC UKB GWAS cups of coffee per day 
#----------------------------------------------exclude people with data in ever using cannabis
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB_CCPD_raw/filePath`; do 
	fileName=`basename $filePath`; 
	echo "fileName=$fileName"; 
	awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $locUKB_CCPD_QC1/${fileName}.allOccurrenceDup; 
done;

# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
## output files are headerless
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB_CCPD_raw/filePath`; do 
	fileName=`basename $filePath`; 
	echo "fileName=$fileName";
	awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $locUKB_CCPD_QC2/${fileName}.allOccurrenceDupRemoved;
done

## Store file paths as a file
realpath $locUKB_CCPD_QC2/*.allOccurrenceDupRemoved > $locUKB_CCPD_QC2/filePath

# Further remove non SNPid (rs number) in the SNP field from files from previous step
# Then remove ambiguous SNPs and insertions, deletions
## The input files have 16 columns. Make sure headers in the awk BEGIN match the number
for filePath in `cat $locUKB_CCPD_QC2/filePath`;do 
	fileName=`basename $filePath`; 
	echo "fileName=$fileName";
	awk 'BEGIN {FS="\t";OFS="\t"; print "SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM",
	"P_BOLT_LMM" } {if( $1 ~/rs/) print $0}' $filePath | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${locUKB_CCPD_QC3}/${fileName}.ambiguousSNPRemoved;
done

# Combine 22 autosome GWAS files as a single file
## Save file paths in a file
realpath ${locUKB_CCPD_QC3}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_all_coffee_cpd.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved > ${locUKB_CCPD_QC3}/filePath
ls
## Vertically combine all 22 files excluding their header rows
cat ${locUKB_CCPD_QC3}/filePath | xargs awk '(FNR>1){print $0}'> ${locUKB_CCPD_QC3}/QCed-GWAS-UKB-CCPD

##Write header back to the file
### Change P_BOLT_LMM_INF to PVALUE for LD score regression to interpret the p value column easily
cd $locUKB_CCPD_QC3;
file1=revised_bolt_imputed_ukb_imp_chr1_v3_all_coffee_cpd.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved;
file2=QCed-GWAS-UKB-CCPD;
cat <(head -1 $file1 | sed 's/P_BOLT_LMM/PVALUE/g') <(cat $file2) > QCed-GWAS-UKB-CCPD_headed


#---------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ QC UKB GWAS estimated standard drinks per week 
#----------------------------------------------exclude people with data in ever using cannabis
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB_ESDPW_raw/filePath`; do 
fileName=`basename $filePath`; 
echo "fileName=$fileName"; 
awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $locUKB_ESDPW_QC1/${fileName}.allOccurrenceDup; 
done;

# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
## output files are headerless
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB_ESDPW_raw/filePath`; do 
fileName=`basename $filePath`; 
echo "fileName=$fileName";
awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $locUKB_ESDPW_QC2/${fileName}.allOccurrenceDupRemoved;
done

## Store file paths as a file
realpath $locUKB_ESDPW_QC2/*.allOccurrenceDupRemoved > $locUKB_ESDPW_QC2/filePath

# Further remove non SNPid (rs number) in the SNP field from files from previous step
# Then remove ambiguous SNPs and insertions, deletions
## The input files have 16 columns. Make sure headers match them
for filePath in `cat $locUKB_ESDPW_QC2/filePath`;do 
fileName=`basename $filePath`; 
echo "fileName=$fileName";
awk 'BEGIN {FS="\t";OFS="\t"; print "SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","P_BOLT_LMM" } {if( $1 ~/rs/) print $0}' $filePath | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${locUKB_ESDPW_QC3}/${fileName}.ambiguousSNPRemoved;
done

# Combine 22 autosome GWAS files as a single file
## Save file paths in a file
realpath ${locUKB_ESDPW_QC3}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_complete_alcohol_unitsweekly.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved > ${locUKB_ESDPW_QC3}/filePath

## Vertically combine all 22 files excluding their header rows
cat ${locUKB_ESDPW_QC3}/filePath | xargs awk '(FNR>1){print $0}'> ${locUKB_ESDPW_QC3}/QCed-GWAS-UKB-ESDPW

##Write header back to the file
cd $locUKB_ESDPW_QC3;
file1=revised_bolt_imputed_ukb_imp_chr1_v3_complete_alcohol_unitsweekly.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved;
file2=QCed-GWAS-UKB-ESDPW;
cat <(head -1 $file1 | sed 's/P_BOLT_LMM/PVALUE/g') <(cat $file2) > ${file2}_headed;

#---------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ QC UKB GWAS 20161 packs years of smoking 
#----------------------------------------------exclude people with data in ever using cannabis
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB20161_raw/filePath`; do 
fileName=`basename $filePath`; 
echo "fileName=$fileName"; 
awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $locUKB20161_QC1/${fileName}.allOccurrenceDup; 
done;

# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
## output files are headerless
field_RSNum=1; # filed position of SNP Rs ID
for filePath in `cat $locUKB20161_raw/filePath`; do 
fileName=`basename $filePath`; 
echo "fileName=$fileName";
awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $locUKB20161_QC2/${fileName}.allOccurrenceDupRemoved;
done

## Store file paths as a file
realpath $locUKB20161_QC2/*.allOccurrenceDupRemoved > $locUKB20161_QC2/filePath

# Further remove non SNPid (rs number) in the SNP field from files from previous step
# Then remove ambiguous SNPs and insertions, deletions
## The input files have 16 columns. Make sure headers match them
for filePath in `cat $locUKB20161_QC2/filePath`;do 
fileName=`basename $filePath`; 
echo "fileName=$fileName";
awk 'BEGIN {FS="\t";OFS="\t"; print "SNP","CHR","BP","GENPOS","ALLELE1","ALLELE0","A1FREQ","INFO","CHISQ_LINREG","P_LINREG","BETA","SE","CHISQ_BOLT_LMM_INF","P_BOLT_LMM_INF","CHISQ_BOLT_LMM","P_BOLT_LMM" } {if( $1 ~/rs/) print $0}' $filePath | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${locUKB20161_QC3}/${fileName}.ambiguousSNPRemoved;
done

# Combine 22 autosome GWAS files as a single file
## Save file paths in a file
realpath ${locUKB20161_QC3}/revised_bolt_imputed_ukb_imp_chr{1..22}_v3_merged_pack_years_20161.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved > ${locUKB20161_QC3}/filePath

## Vertically combine all 22 files excluding their header rows
cat ${locUKB20161_QC3}/filePath | xargs awk '(FNR>1){print $0}'> ${locUKB20161_QC3}/QCed-GWAS-UKB-PYOS

##Write header back to the file
cd $locUKB20161_QC3;
file1=revised_bolt_imputed_ukb_imp_chr10_v3_merged_pack_years_20161.bgen.assoc.allOccurrenceDupRemoved.ambiguousSNPRemoved;
file2=QCed-GWAS-UKB-PYOS;
cat <(head -1 $file1 | sed 's/P_BOLT_LMM/PVALUE/g') <(cat $file2) > ${file2}_headed;

#-----------------------------
#--Import Bash functions to repeat the two 
#-----------------------------
source ${locScripts}/MR_step03-02-01_function_QC-BOLT-LMM-GWAS-files.sh; 

# Calling functions to QC UKB BOLT-LMM GWAS files for caffeine consumed per day
copy_files_get_filePaths "${ref_UKB_caffeine_NA20453_GWAS}" "revised_bolt_imputed_ukb_imp_chr*v3_caffeine.per.day.bgen.assoc" "$locUKB_caffeine_raw"
remove_multiple_occurrences_SNPs 1 "$locUKB_caffeine_raw/filePath" "$locUKB_caffeine_QC1" "$locUKB_caffeine_QC2"
remove_nonSNPid_ambiguous_SNPs_insertions_deletions "$locUKB_caffeine_QC2/filePath" "$locUKB_caffeine_QC3" "QCed-GWAS-UKB-caffeine-consumed-per-day"

#---------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ QC UKB GWAS 20453 ever taken cannabis
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Print lines with more than 1 occurrence of SNP rs number using awk two file processing
field_RSNum=3; # filed position of SNP Rs ID
for filePath in `cat $locUKB20453_raw/filePath`; do 
fileName=`basename $filePath`; 
echo "fileName=$fileName"; 
awk -F"\t" -v fieldPosSNPID=$field_RSNum 'NR==FNR {count[$fieldPosSNPID]++;next} count[$fieldPosSNPID]>1' $filePath $filePath > $locUKB20453_QC1/${fileName}.allOccurrenceDup; 
done;

# Extract lines with only 1 occurrence of SNP ID (i.e. remove all occurrences of duplicated SNPs)
## output files are headerless
field_RSNum=3; # filed position of SNP Rs ID
for filePath in `cat $locUKB20453_raw/filePath`; do 
fileName=`basename $filePath`; 
echo "fileName=$fileName";
awk -F"\t" -v fieldPosSNPID=${field_RSNum} '(FNR!=1){seen[$fieldPosSNPID]++; a[++count]=$0; key[count]=$fieldPosSNPID} END {for (i=1;i<=count;i++) if (seen[key[i]] == 1) print a[i]}' $filePath > $locUKB20453_QC2/${fileName}.allOccurrenceDupRemoved;
done

## Store file paths as a file
realpath $locUKB20453_QC2/*.allOccurrenceDupRemoved > $locUKB20453_QC2/filePath

# Further remove non SNPid (rs number) in the SNP field from files from previous step
# Then remove ambiguous SNPs and insertions, deletions
## The input files have 16 columns. Make sure headers match them
for filePath in `cat $locUKB20453_QC2/filePath`;do 
fileName=`basename $filePath`; 
echo "fileName=$fileName";
awk 'BEGIN {FS="\t";OFS="\t"; print "#CHROM","POS","ID","REF","ALT","TEST","OBS_CT","OR","SE","L95","U95","T_STAT","P" } {if( $3 ~/rs/) print $0}' $filePath | sed 's/C\tG/ambiguous/g;s/G\tC/ambiguous/g;s/T\tA/ambiguous/g;s/A\tT/ambiguous/g' | grep -v ambiguous > ${locUKB20453_QC3}/${fileName}.ambiguousSNPRemoved;
done

# Combine 22 autosome GWAS files as a single file
## Save file paths in a file
realpath ${locUKB20453_QC3}/X20453_0_0_recoded_plink-ukb_imp_chr{1..22}_v3.plink2.output.X20453_0_0_recoded_plink.glm.logistic.allOccurrenceDupRemoved.ambiguousSNPRemoved > ${locUKB20453_QC3}/filePath

## Vertically combine all 22 files excluding their header rows
cat ${locUKB20453_QC3}/filePath | xargs awk '(FNR>1){print $0}'> ${locUKB20453_QC3}/QCed-GWAS-UKB-ever-taken-cannabis

##Write header back to the file
cd $locUKB20453_QC3;
file1=X20453_0_0_recoded_plink-ukb_imp_chr1_v3.plink2.output.X20453_0_0_recoded_plink.glm.logistic.allOccurrenceDupRemoved.ambiguousSNPRemoved ;
file2=QCed-GWAS-UKB-ever-taken-cannabis;
cat <(head -1 $file1 ) <(cat $file2) > ${file2}_headed;

##################################################################################################################
##---------------------------------This is the end of this file-------------------------------------------------##
##################################################################################################################
