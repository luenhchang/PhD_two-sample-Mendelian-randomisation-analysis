#!/bin/bash
## file name: MR_step04-01-02_LD-clumping-criteria.sh
## date created: 20190403
## purpose: Create a file with LD clumping criteria

# Type 	File
#--------------------------------------------------------------------------------------------------------------
# Outpu ${locInput}/LD-SNP-clumping-criteria.txt
#---------------------------------------------------------------------------------------------------------------

## Time 	Change
##--------------------------------------------------------------------------------------------------------------
## 20190403	Exported LD-SNP-clumping-criteria.txt
##--------------------------------------------------------------------------------------------------------------

## Locations of main folders
workingDir="/mnt/lustre/working/lab_nickm/lunC";
locMR="${workingDir}/MR_ICC_GSCAN_201806/data";
locGSCAN="$locMR/noICC_results";
locLDBasedClumping=$locGSCAN/LD-based_SNP_clumping;
locInput=$locLDBasedClumping/input;
mkdir -p $locInput;

# Create a here document with 3 rows
## row1: header
## row2 to last row: each clumping criteria 
cd $locInput;
(cat <<- _EOF_
	r2 LDWindow p1 p2 
	0.01 10000 5e-8 1e-6
	0.01 10000 1e-5 1e-5
	_EOF_
) > ${locInput}/LD-SNP-clumping-criteria.txt

##---------------------------------This is the end of this file-------------------------------------------------##

