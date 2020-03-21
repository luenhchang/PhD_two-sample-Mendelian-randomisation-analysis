## File name: MR_step04-04_history_run-COJO_HPC-Utility.sh
## Modified from: ${homeDir}/scripts/PRS_UKB_201711/PRS_UKB_201711_step00-01_history_runGWAS_HPCUtility.sh
## Date created: 20190302
## Run dependency: /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/documentation-run-GCTA-COJO-analysis.xlsx
## Note: 
## GCTA-COJO website: http://cnsgenomics.com/software/gcta/#COJO
#---------------------------------------------------------------------------------------------------------------------------------------------------------
# --cojo-slct	Perform a stepwise model selection procedure to select independently associated SNPs. Results will be saved in a *.jma file with additional file *.jma.ldr showing the LD correlations between the SNPs.
# --cojo-wind 	Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified.
# --cojo-p 		Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction with the option --cojo-slct.
# --cojo-file 	Input the summary-level statistics from a meta-analysis GWAS
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

## Purpose: (1) step by step guideline for using GUI HPC_Utility.jar to run multi-SNP-based conditional & joint association analysis using GWAS summary data (GCTA-COJO) 
## Time 	Change
##----------------------------------------------------------------------------------------------------
## 20190304	submitted job 7063345, 7063354, 7063355, 7063356, 7063357 
## 20190302 submitted job 6990449,6990450,6990451,6990452(dup),6990453
##----------------------------------------------------------------------------------------------------

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Input $locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-ESDPW_headed
# Input $locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB3456_headed
# Input $locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-PYOS_headed
# Input $locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-CCPD_headed


# Outpu $locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/COJO-QCed-GWAS-UKB-ESDPW_headed.sh
# Outpu $locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/

# Outpu $locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/COJO-QCed-GWAS-UKB3456_headed.sh
# Outpu $locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/QCed-GWAS-UKB3456_headed_clumps_kb_10000_maf_0.01_cojo_p_5e-8.jma.cojo

# Outpu $locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/COJO-QCed-GWAS-UKB-PYOS_headed.sh
# Outpu $locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/QCed-GWAS-UKB-PYOS_headed_clumps_kb_10000_maf_0.01_cojo_p_5e-8.jma.cojo

# Outpu $locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/COJO-QCed-GWAS-UKB-CCPD_headed.sh
# Outpu $locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO/QCed-GWAS-UKB-CCPD_headed_clumps_kb_10000_maf_0.01_cojo_p_5e-8.jma.cojo
#---------------------------------------------------------------------------------------------------------------

UKB_GWAS_BOLTLMM_dir="/mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas/BOLT_LMM/"

homeDir="/mnt/backedup/home/lunC";
inputPhenoLabStuartMa="${homeDir}/reference/data/UKBB_500k/versions/lab_stuartma/pheno";
inputPhenoChang="${homeDir}/data/UKBionbank_phenotype";
locScripts="${homeDir}/scripts/MR_ICC_GSCAN_201806/";

workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";
outputMain="${workingDir}/";

# Step 1: Copy file /mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/HPC-Utility/HPC_Utility.jar and /mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/HPC-Utility/lib to D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR-HPC-Utility-201903 

## launch Windows CMD
## change directory to D:\Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR-HPC-Utility-201903
### (1) type "d:"
### (2) type "cd Now\library_genetics_epidemiology_GWAS_largeFiles\QIMR-HPC-Utility-201903
### (3) type "java -jar HPC_Utility.jar"
### (4) Select cojo as utility; Enter QIMR username, password to establish SFTP connection

# Step 2: run HPC_Utility.jar. On HPC_Utility.jar window
## By default, this software generates files in the input file directory, which you have no permission to write files. So Step 1 copies file to your own directory
## Button           Action
##---------------------------------------------------------------------------------------------------------------------------------------------------
## output dir       Directory of output files 
mkdir -p $locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO $locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO $locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO $locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/GCTA-COJO 

mkdir -p $locMR/ICC-cannabis-ever/GCTA-COJO $locMR/noICC_results/GCTA-COJO

## GWAS file:   Directory of input GWAS summary statistic file
# Conso	Trait	Directory 
#----------------------------------------------------------------------------------------------------------------------------------------------------
# UKB 	ESDPW 	$locMR/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-ESDPW_headed
# UKB	3456CPD	$locMR/UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB3456_headed
# UKB	20161PYOS $locMR/UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-PYOS_headed
# UKB	CCPD	$locMR/UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-CCPD_headed
# GSCAN	SI		$locMR/noICC_results/si_noICC.txt
# ICC	Cannabis use	$locMR/ICC-cannabis-ever/Cannabis_ICC_UKB_small.txt
cp -n $locMR/Cannabis_ICC_UKB_small.txt $locMR/ICC-cannabis-ever/Cannabis_ICC_UKB_small.txt
#----------------------------------------------------------------------------------------------------------------------------------------------------

## bfile /reference/data/UKBB_500k/versions/lab_stuartma/LD_reference/LD_ref_201803_rsOnly as default


## SNP 		variant ID (RS number or CHR:BP?)
## Allele1	effect allele (=alternative allele)
## Allele0	the other allele (i.e. reference allele)
## A1freq	frequency of the effect allele Allele1 (=?MAF)
## MAF		minor allele frequency (i.e. frequency of Allele1 in the population; GSCAN cannot provide MAF)
## BETA		effect size. For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error
## SE		standard error
## P 		p value of association between a phenotype and genotypes

## Parameter 	Value
##----------------------
## kb			10000
## MAF			0.01
## cojo-p		5e-8
## #samples		manually enter overall sample size
## ncpu			8
## mem			65000mb
## time			55 hours
##----------------------

## Hit "generate script" key

## Results *_clumps_kb_10000_maf_0.01_cojo_p_5e-8.jma.cojo
#-----------------------------------
# Conso 	Trait 		Numb-SNPs	SNPs
# ICC+UKB 	CU			4			rs1368740, rs4099556, rs9919557, rs17761723.
# UKB		3456CPD		3			rs112878080,rs56113850, rs67210567
# UKB		20161PYOS	21			rs77244502	rs2867113	rs4953152	rs12996374	rs13393501	rs1561114	rs4894637	rs56369810	rs6942129	rs10274594	rs215670	rs6962772	rs1419748	rs373682067	rs17309874	rs112282219	rs577011502	rs8042849	rs12599637	rs2316205	rs59145036
# UKB		CCPD		21
# UKB		ESDPW		39
#-----------------------------------
 


##----------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------#
#--------------------------Document history of running GWAS using HPC Utility
#------------------------------------------------------------------------------------------------------------------------------#

# Parameter			Value
#-------------------------------------------------------------------------------------------------------------------------------
# tool				BOLT-LMM
# output dir        /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3456_numCigareDaily
# pheno file        /mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb3456_numCigareDaily/ukb3456.phenoUtility
# pheno 			3456.0.0
# exclude sample    /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
# qsub jobs			4827782..4827803 (These jobs took > 21:30 : 4827788-4827794)
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3456_numCigareDaily/BOLT-LMM-ukb3456.phenoUtility_output_3456-0.0/revised_bolt_imputed_HRC.chr{1..22}.bgen.assoc
cp -n ${outputMain}/GWAS_UKB3456_numCigareDaily/ukb3456.phenoUtility-BOLT-LMM_3456-0.0.sh ${locScripts}/PRS_UKB_201711_step00-02_run_GWAS_BOLT-LMM_UKB3456_numCigareDaily.sh
 
# tool				plink2
# output dir		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20160_everSmoked
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb20160_everSmoked/ukb20160.phenoUtility.recoded
# pheno				X20160_recode
# conti/binary		binary
# exclude samples   as default 
# qsub jobs			4827689..4827710
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20160_everSmoked/plink2-ukb20160.phenoUtility.recoded_output_X20160_recode/X20160_recode_HRC.chr{1..22}.bgen.X20160_recode.glm.logistic
cp -n ${outputMain}/GWAS_UKB20160_everSmoked/smoking.pheno-plink2_ever_smoked.sh ${locScripts}/PRS_UKB_201711_step00-03_run_GWAS_plink2_UKB20160_everSmoked.sh

# tool				plink2
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb20116_smokingStatus/ukb20116.phenoUtility.recoded
# pheno				X20116_recodeFinal
# output dir		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20116_smokingStatus
# conti/binary		binary
# qsub jobs			4828006..4828027
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB20116_smokingStatus/plink2-ukb20116.phenoUtility.recoded_output_X20116_recodeFinal/X20116_recodeFinal_HRC.chr{1..22}.bgen.X20116_recodeFinal.glm.logistic

# tool				BOLT-LMM
# output dir 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3436_ageStartedSmokingInCurrentSmokers
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb3436_ageStartedSmokingInCurrentSmokers/ukb3436.phenoUtility.recoded
# pheno 			X3436_recodeMean
# exclude sample    /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
# qsub jobs			4828028..4828049
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB3436_ageStartedSmokingInCurrentSmokers/BOLT-LMM-ukb3436.phenoUtility.recoded_output_X3436_recodeMean/revised_bolt_imputed_HRC.chr{1..22}_X3436_recodeMean.bgen.assoc

# tool				BOLT-LMM
# output dir 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB1559_numberStandardDrinksPerWeek
# pheno file		/mnt/backedup/home/lunC/data/UKBionbank_phenotype/ukb1559_numberStDrinksPerWeek/alcohol.recoded.weeklyunits.full.pheno

# pheno 			NSDPW
# exclude sample    /reference/data/UKBB_500k/versions/lab_stuartma/exclude_samples/430k_whiteBrit_default.exclude.list
# qsub jobs			4828050..4828071
# GWAS files 		/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB1559_numberStandardDrinksPerWeek/BOLT-LMM-alcohol.recoded.weeklyunits.full.pheno_output_NSDPW/revised_bolt_imputed_HRC.chr{1..22}_NSDPW.bgen.assoc
#--------------------------------------------------------------------------------------------------------------------------------

# Create output folders
mkdir -p $outputMain/GWAS_UKB20160_everSmoked $outputMain/GWAS_UKB2907_everStoppedSmokingFor6+Months $outputMain/GWAS_UKB3436_ageStartedSmokingInCurrentSmokers $outputMain/GWAS_UKB1559_numberStandardDrinksPerWeek;  

#---------------------------------------------------------------------------------------------------------------#
#------------------------------copy output GWAS file directory to a CSV file
#---------------------------------------------------------------------------------------------------------------#

# file path: /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_info.csv

# Copy this file for similar job
cp -n ${homeDir}/scripts/MR_ICC_GSCAN_201806/MR_step04-04_history_run-COJO_HPC-Utility.sh ${homeDir}/scripts/MR_ICC_GSCAN_201806/MR_step08-01_history_run-LDSC-genetic-correlation_HPC-Utility.sh
##---------------------------------This is the end of this file-------------------------------------------------##