##########################################################################################################
# filename: MR_step05-02_calculate-proportion-of-trait-variation-explained-by-genetic-instruments.R
# program author: Chang
# purpose: Calculate proportion of exposure variance explained by clumped SNPs (R2). Compare number of exposure SNPs (independent ones), total R2 with GSCAN studies and other studies
# date created: 20191024
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Input paste0(loc.clumped.GSCAN,"GWAS_from-clumped-SNPs_si_noICC_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
# Input /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/ICC-cannabis-ever/QC4_GWAS_from_clumped_SNPs/GWAS_from-clumped-SNPs_GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6
# Input /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6_sample-size-added
# Input /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6_sample-size-added
# Input /mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB_imputed201803/UKB20116_smokingStatus/QC4_subsetColumns/GWAS_UKB_SS
# Input /mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB20453-ever-taken-cannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-ever-taken-cannabis_headed
# Input /mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20160_everSmoked/ukb20160.phenoUtility.recoded
# Input /mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded
# Input /reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas/BOLT_LMM/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/phenotype
# Input /reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas/BOLT_LMM/UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/phenotype

# Outpu paste0(loc.instrument.strength,"proportion-MR-exposure-variance-explained-by-clumped-SNPs.tsv")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20191024  Exported paste0(loc.instrument.strength,"proportion-MR-exposure-variance-explained-by-clumped-SNPs.tsv")
#-----------------------------------------------------------------------------------------

#-------------------------------------------
# Folder locations under Stuart MC lab
#-------------------------------------------
lab_stuartmaDir <- "/reference/data/UKBB_500k/versions/lab_stuartma/"
loc.BOLTLMM <- paste0(lab_stuartmaDir,"draft_gwas/BOLT_LMM/") 

#-------------------------------------------
# Folder locations under my home directory
#-------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locPlots_MR <- paste0(homeDir,"plots/MR_ICC_GSCAN_201806/")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))

#--------------------------------------------------------------------
# Folder locations under my working 
#--------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locICC <- paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
loc.clumped.GSCAN <- paste0(locICC,"noICC_results/QC4_GWAS_from_clumped_SNPs/")
loc.instrument.strength <- paste0(workingDir,"MR_ICC_GSCAN_201806/instrument-strength/")

#--------------------------------------------------------------------
# Import clumped exposure GWASs
#--------------------------------------------------------------------
# Clumped SNPs for GSCAN exposure regular smoking initiation
GSCAN.SI.5e8.small <- ImportASpaceSeparatedFile(input.file.path = paste0(loc.clumped.GSCAN,"GWAS_from-clumped-SNPs_si_noICC_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
                          ,data.name = "GSCAN.SI.5e8") %>% 
  dplyr::select_(.dots=c("SNP","REF","ALT","BETA","PVALUE")) # dim(GSCAN.SI.5e8) 8 5

ICC.CI.5e8.small <- ImportASpaceSeparatedFile(input.file.path = "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/ICC-cannabis-ever/QC4_GWAS_from_clumped_SNPs/GWAS_from-clumped-SNPs_GWAS-ICC-CI_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6"
                                              ,data.name = "ICC.CI.5e8") %>%
  dplyr::select_(.dots=c("SNP","Allele1","Allele2","MAF","Effect","P")) # dim(ICC.CI.5e8.small) 4 6

UKB.caffeine.5e8.small <- ImportASpaceSeparatedFile(input.file.path = "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/GWAS_from-clumped-SNPs_GWAS-UKB-caffeine_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6_sample-size-added", data.name = "UKB.caffeine.5e8.small") %>% 
  dplyr::select_(.dots=c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA","PVALUE")) # dim(UKB.caffeine.5e8.small) 22 6

UKB.ESDPW.5e8.small <- ImportASpaceSeparatedFile(input.file.path = "/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/GWAS_from-clumped-SNPs_GWAS-UKB-ESDPW_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6_sample-size-added",data.name = "UKB.ESDPW.5e8") %>%
  dplyr::select_(.dots=c("SNP","ALLELE1","ALLELE0","A1FREQ","BETA","PVALUE")) # dim(UKB.ESDPW.5e8.small) 39 6

#----------------------------------------
# Import UKB GWAS for BMI
## This file is to get effect allele frequency to use as MAF in R2 calculation if effect allele frequency is unavaiable in GSCAN
#----------------------------------------
# GWAS for UKB BIM 
ImportATabSeparatedFile(input.file.path = paste0(loc.BOLTLMM,"BMI/BMI.assoc")
                        ,data.name = "GWAS.UKB.BMI") # dim(GWAS.UKB.BMI) 10272963       16

GWAS.UKB.BMI <- GWAS.UKB.BMI %>% 
  dplyr::select_(.dots=c("SNP","ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM")) # dim(GWAS.UKB.BMI) 10272963 5

# GWAS for 20160 ever smoked
GWAS.UKB.20160 <- ImportATabSeparatedFile(input.file.path = "/mnt/lustre/working/lab_nickm/lunC/PRS_UKB_201711/GWASSummaryStatistics/GWAS_UKB_imputed201803/UKB20116_smokingStatus/QC4_subsetColumns/GWAS_UKB_SS"
                                          ,data.name = "GWAS.UKB.20160") %>%
  dplyr::select_(.dots=c("SNPID","BETA")) %>%
  dplyr::rename(SNP=SNPID) # dim(GWAS.UKB.20160) 5093954       2

# Read a TSV file with # in header row
GWAS.UKB.20453 <- read.table(file="/mnt/lustre/working/lab_nickm/lunC/MR_ICC_GSCAN_201806/data/UKB20453-ever-taken-cannabis/QC3_remove_ambiguousSNPs_indel/QCed-GWAS-UKB-ever-taken-cannabis_headed", header=TRUE, comment.char = "", sep="\t", check.names = FALSE, stringsAsFactors = FALSE) %>%
  dplyr::mutate(BETA= log(OR)) %>%
  dplyr::select_(.dots=c("ID","OR","BETA")) %>%
  dplyr::rename(SNP=ID) # dim(GWAS.UKB.20453) 6021166       3

#--------------------------------------------------------------------------
# Import phenotype files to calculate variance of phenotypes
#--------------------------------------------------------------------------
UKB.pheno.20160.small <- ImportASpaceSeparatedFile(input.file.path = "/mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20160_everSmoked/ukb20160.phenoUtility.recoded"
                          ,data.name = "UKB.pheno.20160") %>%
  dplyr::select_(.dots=c("FID","IID","X20160_recode")) # dim(UKB.pheno.20160.small) 487409  3

UKB.pheno.20453.small <- ImportASpaceSeparatedFile(input.file.path = "/mnt/backedup/home/lunC/data/UKBiobank_phenotype/ukb20453_everTakenCannabis/ukb20453.phenoUtility.recoded"
                          ,data.name = "UKB.pheno.20453") %>% 
  dplyr::select_(.dots=c("FID","IID","X20453_0_0_recoded_plink")) # dim(UKB.pheno.20453.small) 487184 3

UKB.pheno.caffeine.small <- ImportASpaceSeparatedFile(input.file.path = "/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas//BOLT_LMM/UKB-estimated-caffeine-consumed-per-day-thru-regular-coffee-and-tea_IID-NA-in-UKB20453-everUsedCannabis/phenotype", data.name = "UKB.pheno.caffeine") %>%
  dplyr::select_(.dots = c("FID","IID","caffeine.per.day")) # dim(UKB.pheno.caffeine.small) 333683      3

UKB.pheno.ESDPW.small <- ImportASpaceSeparatedFile(input.file.path = "/reference/data/UKBB_500k/versions/lab_stuartma/draft_gwas//BOLT_LMM/UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/phenotype", data.name = "UKB.pheno.ESDPW") %>%
  dplyr::select_(.dots=c("FID","IID","complete_alcohol_unitsweekly")) # dim(UKB.pheno.ESDPW.small) 333683      3
#--------------------------------------------------------------------------
# Calculate the phenotypic variation
#--------------------------------------------------------------------------
varia.UKB.pheno.20160 <- var(UKB.pheno.20160.small$X20160_recode, na.rm = TRUE) # 0.2403426
varia.UKB.pheno.20453 <- var(UKB.pheno.20453.small$X20453_0_0_recoded_plink, na.rm = TRUE) # 0.1717352
varia.UKB.pheno.caffeine <- var(UKB.pheno.caffeine.small$caffeine.per.day, na.rm = TRUE) # 28603.29
varia.UKB.pheno.ESDPW <- var(UKB.pheno.ESDPW.small$complete_alcohol_unitsweekly, na.rm = TRUE) # 270.6388

#--------------------------------------------------------------------------
# Calculate R2 by individual SNPs and total R2
#--------------------------------------------------------------------------
# Left join clumped SNP RSID, effect allele freq, and beta
t1 <- list( GSCAN.SI.5e8.small[,c("SNP","PVALUE")]
           ,GWAS.UKB.BMI[,c("SNP","ALLELE1","ALLELE0","A1FREQ")]
           ,GWAS.UKB.20160) %>% 
  purrr::reduce(dplyr::left_join, by=c("SNP")) %>%
  dplyr::mutate(R.sqaured.SNP= dplyr::case_when( BETA==NA~ as.numeric(NA)
                                                ,BETA== -Inf ~ as.numeric(NA)
                                                ,TRUE ~ 2*A1FREQ*(1-A1FREQ)*(BETA**2)/varia.UKB.pheno.20160
                                                )) # dim(t1) 8 7
# Append the summed R-squred 
t1.total <- data.frame( SNP="Total"
                       ,PVALUE=NA
                       ,BETA=NA
                       ,ALLELE1="NA"
                       ,ALLELE0="NA"
                       ,A1FREQ=NA
                       ,R.sqaured.SNP=sum(t1$R.sqaured.SNP, na.rm = TRUE))

t1 <- rbind(t1, t1.total) # dim(t1) 9 7
#            SNP   PVALUE   A1FREQ        BETA R.sqaured.SNP
# 1    rs2155646 1.03e-10 0.612520          NA            NA
# 2    rs7613360 1.88e-09 0.602016 -0.01839410  0.0006745746
# 3   rs62025923 2.60e-09 0.795548          NA            NA
# 4    rs2162965 2.79e-09 0.751556 -0.03636630  0.0020548843
# 5    rs6756212 4.73e-09 0.461532  0.00740253  0.0001133239
# 6  rs184584210 1.06e-08       NA          NA            NA
# 7     rs325535 2.02e-08 0.837055        -Inf            NA
# 8     rs951740 4.58e-08 0.377472 -0.02853530  0.0015922373
# 9        Total       NA       NA          NA            NA
# 10       Total       NA       NA          NA  0.0044350201

# ICC cannabis initiation 
t2 <- list(ICC.CI.5e8.small[,c("SNP","P")]
          ,GWAS.UKB.BMI[,c("SNP","ALLELE1","ALLELE0","A1FREQ")]
          ,GWAS.UKB.20453[,c("SNP","BETA")]) %>%
  purrr::reduce(dplyr::left_join, by=c("SNP")) %>% 
  dplyr::mutate(R.sqaured.SNP= 2*A1FREQ*(1-A1FREQ)*(BETA**2)/varia.UKB.pheno.20453) %>%
  dplyr::rename(PVALUE=P) # dim(t2) 4 7

# Append the summed R-squred 
t2.total <- data.frame(SNP="Total",PVALUE=NA,ALLELE1="NA",ALLELE0="NA",A1FREQ=NA,BETA=NA
                       ,R.sqaured.SNP=sum(t2$R.sqaured.SNP, na.rm = TRUE))
t2 <- rbind(t2, t2.total)

# This used BETA from the UKB GWAS
#          SNP Allele1 Allele2    MAF         P        BETA R.sqaured.SNP
# 1  rs1368740       a       g 0.7510 1.338e-13 -0.06728255   0.009858588
# 2  rs9919557       t       c 0.6134 6.534e-11          NA            NA
# 3  rs4099556       a       g 0.8242 5.908e-09 -0.07152393   0.008632265
# 4 rs17761723       t       c 0.3453 2.153e-08 -0.03610605   0.003432175
# 5      Total      NA      NA     NA        NA          NA   0.021923029

# This used Effect from the ICC GWAS
#          SNP Allele1 Allele2    MAF  Effect         P R.sqaured.SNP
# 1  rs1368740       a       g 0.7510 -0.0768 1.338e-13   0.012844946
# 2  rs9919557       t       c 0.6134 -0.0610 6.534e-11   0.010276278
# 3  rs4099556       a       g 0.8242  0.0699 5.908e-09   0.008244729
# 4 rs17761723       t       c 0.3453  0.0528 2.153e-08   0.007339684
# 5      Total      NA      NA     NA      NA        NA   0.038705638

# UKB caffeine consumption
t3 <- UKB.caffeine.5e8.small %>%
  dplyr::mutate(R.sqaured.SNP= 2*A1FREQ*(1-A1FREQ)*(BETA**2)/varia.UKB.pheno.caffeine) %>%
  dplyr::select_(.dots=c("SNP","PVALUE","ALLELE1","ALLELE0","A1FREQ","BETA","R.sqaured.SNP"))

t3.toal <- data.frame(SNP="Total",PVALUE=NA,ALLELE1="NA", ALLELE0="NA",A1FREQ=NA,BETA=NA
                      ,R.sqaured.SNP=sum(t3$R.sqaured.SNP, na.rm = TRUE))
t3 <- rbind(t3, t3.toal) # dim(t3) 24 7


# UKB ESDPW
t4 <- UKB.ESDPW.5e8.small %>%
  dplyr::mutate(R.sqaured.SNP= 2*A1FREQ*(1-A1FREQ)*(BETA**2)/varia.UKB.pheno.ESDPW) %>%
  dplyr::select_(.dots=c("SNP","PVALUE","ALLELE1","ALLELE0","A1FREQ","BETA","R.sqaured.SNP")) # dim(t4) 39 7

t4.total <- data.frame(SNP="Total",PVALUE=NA,ALLELE1="NA", ALLELE0="NA",A1FREQ=NA,BETA=NA
                       ,R.sqaured.SNP=sum(t4$R.sqaured.SNP, na.rm = TRUE))

t4 <- rbind(t4, t4.total) # dim(t4) 40 7

# rbind all the data sets
t1$exposure.trait <- "SI"
t2$exposure.trait <- "CI"
t3$exposure.trait <- "caffeine"
t4$exposure.trait <- "ESDPW"

all.t <- dplyr::bind_rows(t1,t2,t3,t4) %>%
  dplyr::select_(.dots=c("exposure.trait","SNP","PVALUE","ALLELE1","ALLELE0","A1FREQ","BETA","R.sqaured.SNP"))

ExportFileTabSeparated(data=all.t
                       , missing.values.as = ""
                       , output.file.path = paste0(loc.instrument.strength,"proportion-MR-exposure-variance-explained-by-clumped-SNPs.tsv"))

# Compare number of indepdent SNPs and total R2 with GSCAN's number of indepdent SNPs and heritability




