# ------------------------------------------------------------------------------------------------
# Program       : MR_step03-03_compare-sample-sizes_QC-GWAS-GSCAN_different-samples.R
# Modified from : 
# Date created  : 20181121
# Purpose       : List clumped SNPs in GSCAN noICC and GSCAN previous sample
# Note          : 
#----------------------------------------------------------------------------------------------
# Run dependency: 
# Type File
#---------------------------------------------------------------------------------------------
# Input paste0(,"")
# Outpu paste0(,"")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 201811
#----------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

dir.working <- "/mnt/lustre/working/lab_nickm/lunC/"
loc.GSCAN <- paste0(dir.working,"MR_ICC_GSCAN_201806/data/")
#loc.GSCAN <- "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/"
#loc.no23andMe <- paste0(loc.GSCAN,"no23andMe_results/QC3_remove_ambiguousSNPs_indel/")
loc.clumped.no23andMe <- paste0(loc.GSCAN,"no23andMe_results/LD-based_SNP_clumping/output/")
#loc.noICC <- paste0(loc.GSCAN,"noICC_results/QC3_remove_ambiguousSNPs_indel/")
loc.clumped.noICC <- paste0(loc.GSCAN,"noICC_results/LD-based_SNP_clumping/output/")

# Create a two column file with directory to clumped GWAS files
#file.path.GWAS.no23andMe <- Sys.glob(paste0(loc.no23andMe,"*_no23andme.ambiguousSNPRemoved")) # 5 files
file.path.clumped.GWAS.no23andMe <- Sys.glob(paste0(loc.clumped.no23andMe,"*_no23andme_LDWindow-kb-*_R2-0.01_p1-1e-6_p2-1e-6.clumped")) # 10 files 

#file.path.GWAS.noICC <- Sys.glob(paste0(loc.noICC,"*_noICC.ambiguousSNPRemoved")) # 5 files
file.path.clumped.GWAS.noICC <- Sys.glob(paste0(loc.clumped.noICC,"*_LDWindow-kb-*_R2-0.01_p1-1e-6_p2-1e-6.clumped")) # 10 files

path.output <- "/mnt/lustre/reference/data/UKBB_500k/versions/lab_stuartma/collab/lunC/MR_ICC_GSCAN_201806/data/clumpedSNPs.compare.no23andMe.noICC/"

dir.create(path.output)

# Full join no23andMe and noICC
#common.columns <- c("RSID","N","EFFECTIVE_N","Number_of_Studies")
common.columns <- c("CHR","F","SNP","BP","P")

for (i in 1:length(file.path.clumped.GWAS.no23andMe)){
  file.path.no23andMe <- file.path.clumped.GWAS.no23andMe[i]
  file.path.noICC <- file.path.clumped.GWAS.noICC[i]
  
  file.no23andMe <- read.table(file.path.no23andMe
                               , header=T
                               ,sep = ""
                               ,stringsAsFactors = F) %>% 
                      select_(.dots=common.columns)
  
  file.noICC <- read.table(file.path.noICC,header=T
                           ,sep = ""
                           ,stringsAsFactors = F) %>% 
                  select_(.dots=common.columns)
  # No SNPs common to both samples
  #file.left.join <- dplyr::left_join(file.no23andMe,file.noICC, by=c("SNP"="SNP"))
  
}  






