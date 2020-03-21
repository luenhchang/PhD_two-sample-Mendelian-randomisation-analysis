#!/usr/bin/Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step00-02_calculate_cases_controls.R
# Modified from : 
# Date created  : 20180823
# Purpose       : get number of cases and controls, sample prevalence, population prevalence for binary traits from literature
# Run dependency: D:\googleDrive\GoogleDriveAsMaster\Chang_manuscript04_MR-cannabis-licitSubstanceUse\international-journal-of-epidemiology\manuscript4-main-text.docx
# Note: 
#----------------------------------------------------------------------------------------
# Run dependency: 
# Function external: 

# Type  Files
#----------------------------------------------------------------------------------------------
# Input 
# Outpu 
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20190312    Updated calculation of cases and controls using information from manuscript4-main-text.docx
#----------------------------------------------------------------------------------------

# Restrict SNPs to MAF between 0.01 and 0.99
cannabis <-read.table(file="/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/Cannabis_ICC_UKB_small.txt"
                      ,sep=" ",header=T)

summary(cannabis)

subset(cannabis, cannabis$MAF<0.99 & cannabis$MAF>0.01) -> c1

write.table(c1
            ,file="/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/Cannabis_ICC_UKB_small_MAF0.99-0.01"
            ,sep=" "
            ,row.names = F
            , col.names = TRUE
            , quote = F)

#---------------------------------------------------------------------------------------------------------------------
# Calculate number of cases and controls per cannabis initiation (CI) from ICC where UKB is the largest contributing study
## Text below is copied from Methods section of manuscript4-main-text.docx
## Cannabis initiation was assessed by the ICC consortium using self-reported question “Have you taken CANNABIS (marijuana, grass, hash, ganja, blow, draw, skunk, weed, spliff, dope), even if it was a long time ago?” Meta-analysis GWAS summary statistics of this binary phenotype was conducted in the full ICC sample excluding 23andMe (N=164742; 26.76% cannabis users; %females=?; ICC subjects: 35297; UKB subjects: 126785)
#----------------------------------------------------------------------------------------
overall.sample.size.CI <- 164742
prevalence.CI <- 0.2676
numb.cases.CI <- overall.sample.size.CI*prevalence.CI # 44084.96
numb.ctrl.CI <- overall.sample.size.CI*(1-prevalence.CI) # 120657

#---------------------------------------------------------------------------------------------------------------------
# Calculate number of cases and controls per binary trait from GSCAN
## Get sample size information from personal communication with Mengzhen Liu. Text below is copied from Methods section of manuscript4-main-text.docx :
## smoking initiation (SI; N=207 726; 54% smokers; 65.9% females); age of initiation of regular smoking (AI; N=119 239, 58.4% females); cigarettes per day (CPD; N=122 027; 61.5%females); smoking cessation (SC; N=125362; 34% current smokers; 51.9% females); and drinks per week (DPW; N=185 828; 55.15%females).
#----------------------------------------------------------------------------------------
# Copy overall sample size (as contrast to the sample size varying with SNPs)
overall.sample.size.GSCAN.SC <- 125362 
overall.sample.size.GSCAN.SI <- 207726 

# Calculate number of cases in each binary trait using provided prevalence
numb.cases.GSCAN.SC <- overall.sample.size.GSCAN.SC*(1-0.34) # 82738.92
numb.ctrl.GSCAN.SC <- overall.sample.size.GSCAN.SC-numb.cases.GSCAN.SC # 42623.08

numb.cases.GSCAN.SI <- overall.sample.size.GSCAN.SI*0.54 # 112172
numb.ctrl.GSCAN.SI <- overall.sample.size.GSCAN.SI-numb.cases.GSCAN.SI # 95553.96

#--------------------------------------------------------------------------------------------------
# Get population prevalence of binary traits
## Ever stopped smoking (Smoking Cessation) : calculated from http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=2907
num_yes=55437
num_no=72997
mean_prevalence_SC=num_yes/(num_yes+num_no)

#----------------------------------------------------------------


