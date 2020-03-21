#!/bin/bash
## file name: MR_step04-01-01_jobScript_LD-based-SNP-clumping.sh
## old file name: 
## modified from: 
## date created: 20180719
## purpose: select independent SNPs
## Run dependency: 
## How to run this script: 
## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20180719	qsub jobs 
##--------------------------------------------------------------------------------------------------------------

# Pass qsub variables to shell variables
bfile=$v_bfilePath;
filePath=$v_clumpFilePath;
clumping_r2_threshold=$v_clumping_r2_threshold;
clumping_LD_window=$v_clumping_LD_window;
p_threshold_1=$v_p_threshold_1;
p_threshold_2=$v_p_threshold_2;
outputDir=$v_outputDir;
outputResultFileName=$v_outputResultFileName;
plink=$v_plink;

homeDir="/mnt/backedup/home/lunC";
locScripts=${homeDir}/scripts/MR_ICC_GSCAN_201806

#-------------------------------------------------------------------------------------------#
#------------------------------LD-based clumping using plink--------------------------------#
#-------------------------------------------------------------------------------------------#
# The --clump command is used to specify one or more result files (i.e. precomputed analyses of some kind). By default, PLINK scans these files and extracts fields with the headers "SNP" and "P"

module load plink/1.90b4.1

## Clump SNPs 
plink	--bfile ${bfile} \
		--clump ${filePath} \
        --clump-p1 ${p_threshold_1} \
        --clump-p2 ${p_threshold_2} \
        --clump-r2 ${clumping_r2_threshold} \
        --clump-kb ${clumping_LD_window} \
        --out ${outputDir}/${outputResultFileName}
#cp -n ${locScripts}/MR_step04-01-01_jobScript_LD-based-SNP-clumping.sh ${locScripts}/MR_step08-02-02_jobScript_calculate-genetic-correlation-LDSC.sh  		
##---------------------------------This is the end of this file-------------------------------------------------##