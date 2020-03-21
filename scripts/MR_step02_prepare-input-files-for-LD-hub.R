##################################################################################################
# program name      : MR_step02_prepare-input-files-for-LD-hub.R
# modifiied from    : 
# purpose           : 
# programmer  	    : Chang
# date created	    : 20180713
# external function : nil
# Internal function : 
# note			    : 
#---------------------------------------------------------------------------------------
# run dependency  : 
#                   
# Type  File
#---------------------------------------------------------------------------------------
# Input ${locLD}/w_hm3.noMHC.snplist
# Input ${locICC}/Cannabis_ICC_UKB.txt
# Outpu ${locICC}/Cannabis-ICC-UKB_merged-LD-hub-SNP-list.zip
#------------------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20180424  Exported the 2 file above
#----------------------------------------------------------------------------------------
# Locations of main folders
homeDir="/mnt/backedup/home/lunC/";
workingDir="/mnt/lustre/working/lab_nickm/lunC/";

# Folders under the main folders
locICC=paste0(workingDir,"MR_ICC_GSCAN_201806/data/")
locLD=paste0(homeDir,"LD-hub/")

#-------------------------------------------------------------------------------------#
#--------Import LD-hub SNP list and ICC GWAS on cannabis use--------------------------#
#-------------------------------------------------------------------------------------#
# Merge LD hub SNP list with cannabis GWAS
file1=read.table(paste0(locLD,"w_hm3.noMHC.snplist"),header = T, sep="\t", stringsAsFactors = F) # 1215001 obs. of  3 variables
file2=read.table(paste0(locICC,"Cannabis_ICC_UKB.txt"),header = T, sep=" ",stringsAsFactors = F) # 11535592 obs. of  11 variables:

# Left-join file1 (left table) and file2 (right table) with merging key as same column names
## keep columns of right table to those wanted
file2_small=file2[,c("SNP","Allele1","Allele2","MAF","Effect","P","N")]

leftJoin=merge(file1,file2_small,by="SNP",all.x=TRUE) # 1215001 obs. of  9 variables

# Filter out SNPs with missing values in Allele1
variables_to_select=c("SNP","Allele1","Allele2","MAF","Effect","P","N")
leftJoin2=subset(leftJoin,(!is.na(leftJoin[,"Allele1"])),select=variables_to_select) # 1203650 obs. of  7 variables

# Add number of controls
N=164742
leftJoin2$N_ctrl=with(leftJoin2,164742-N)

# Rename columns that are consistent with an example log file
##  Old-colname  New-colname  Def
##---------------------------------------------------------------
##  SNP           SNP         Variant ID (e.g., rs number)
##  Allele1       EA          Allele 1 interpreted as ref allele for signed sumstat.
##  Allele2       NEA         Allele 2, interpreted as non-ref allele for signed sumstat.
##  MAF           EAF         Effective allelic frequency
##  Effect        BETA      beta, OR or Z-score. Must have a direction (i.e. positve or negative sign)
##  P             P           p-Value
##  N_ctrl        N_ctrl      Number of controls, calculated as 164742- N (case numbers). 164742 is the sum of substudy samples in samples.xlsx
##---------------------------------------------------------------

leftJoin3=leftJoin2[,-7]

colnames(leftJoin3)=c("SNP","EA","NEA","EAF","BETA","P-value","N")

# Save output file as a txt file and then zip the txt file
outFilePath=paste0(locICC,"Cannabis-ICC-UKB_merged-LD-hub-SNP-list.txt")
write.table(leftJoin3 # object name the file to export
            ,col.names=T   # keep column names
            ,row.names = F # remove row number
            ,file=outFilePath
            ,dec="."
            ,sep=" "
            ,quote=FALSE
            ,na = "NA" ) # mark missing values as NA

## Zip the txt file under the same folder
setwd(locICC)
zip(zipfile = 'Cannabis-ICC-UKB_merged-LD-hub-SNP-list', files = 'Cannabis-ICC-UKB_merged-LD-hub-SNP-list.txt')

# Next step: 
## copy file to D:\Now\library_genetics_epidemiology_GWAS_largeFiles\international_cannabis_consortium\data
## upload the zip file to http://ldsc.broadinstitute.org/upload_file/#
