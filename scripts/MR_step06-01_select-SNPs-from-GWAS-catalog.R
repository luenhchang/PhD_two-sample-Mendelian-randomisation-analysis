##################################################################################
# filename: MR_step06-01_select-SNPs-from-GWAS-catalog.R
# program author: Chang
# purpose: is there observational association between exposure traits and outcome trait, giving support for doing MR
# date created: 20180825
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Input paste0(locICCRef,"top-hits-substance-use.txt")

# Outpu published_topSNPs_alcohol
# Outpu published_topSNPs_coffee
# Outpu published_topSNPs_nicotine
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20180
#-----------------------------------------------------------------------------------------

# Folder locations under my home directory
homeDir="/mnt/backedup/home/lunC/";

# Folder locations under Stuartma lab
stuartmaDir="/reference/data/UKBB_500k/versions/lab_stuartma/"
stuartma_pheno=paste0(stuartmaDir,"pheno/")
ref_UKB_GWAS_BOLT_LMM=paste0(stuartmaDir,"draft_gwas/BOLT_LMM/")
ref_UKB_ESDPW=paste0(ref_UKB_GWAS_BOLT_LMM,"UKB_estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/")
ref_UKB_CCPD=paste0(ref_UKB_GWAS_BOLT_LMM,"UKB-cups-of-coffee-per-day_IID-NA-in-UKB20453-everUsedCannabis/")
ref_UKB20161=paste0(ref_UKB_GWAS_BOLT_LMM,"UKB20161-pack-years-of-smoking_IID-NA-in-UKB20453-everUsedCannabis/")

# Folders under lunC home  
locScripts=paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
locUKB3456_pheno=paste0(homeDir,"data/UKBionbank_phenotype/ukb3456_numCigareDaily/")

# Folders under lunC working
workingDir="/mnt/lustre/working/lab_nickm/lunC/";
locICCRef=paste0(workingDir,"MR_ICC_GSCAN_201806/reference/")
loc.GWAS.catalog <- paste0(workingDir,"MR_ICC_GSCAN_201806/GWAS-catalog-download/")

locICC=paste0(workingDir,"MR_ICC_GSCAN_201806/data/") # location of outcome data
#locGSCAN=paste0(locICC,"no23andMe_results/QC4_GWAS_from_clumped_SNPs/LDWindow10000kb/") # location of exposure data files
locGSCAN.QC3 <- paste0(locICC,"noICC_results/QC3_remove_ambiguousSNPs_indel/")
locGSCAN.QC4 <- paste0(locICC,"noICC_results/QC4_GWAS_from_clumped_SNPs/")

locUKB3456_QC4=paste0(locICC,"UKB3456-numCigareDaily_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")
#locUKB_ESDPW_QC4=paste0(locICC,"UKB-estimated-standard-drinks-per-week/QC4_GWAS_from_clumped_SNPs/")
locUKB_ESDPW_QC4=paste0(locICC,"UKB-estimated-standard-drinks-per-week_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")
locUKB_CCPD_QC4=paste0(locICC,"UKB-cups-coffee-per-day_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")
locUKB20161_QC4=paste0(locICC,"UKB20161-packs-years-of-smoking_IID-NA-in-UKB204534-everUsedCannabis/QC4_GWAS_from_clumped_SNPs/")

locPlots_MR=paste0(homeDir,"plots/MR_ICC_GSCAN_201806/")

library(TwoSampleMR)
library(dplyr)
library(stringr)


#--------------------------------------------------------------------------------------------------------
# Find SNPs with proven or plausible biological effects on the target exposure
## SNPs selected from Tables of Prom-Wormley et al 2017 The genetic epidemiology of substance use disorder A review.pdf
## with p values <5e-08
#--------------------------------------------------------------------------------------------------------
# SNPs selected from Prom-Wormley et al 2017's tables nrow(topHits_drug)
cls <- c(pValue="numeric")

# Warning: this input file top-hits-substance-use.txt doesn't have six comma separated values in every line
# topHits_drug= read.table(paste0(locICCRef,"top-hits-substance-use.txt")
#                          ,sep=","
#                          ,header=T
#                          ,stringsAsFactors = F
#                          ,colClasses=cls) %>%
#   filter(pValue < 5e-08) %>%
#   filter(! EffectSize %in% c("Not Reported","?")) %>%
#   filter(EffectAllele %in% c("A","T","C","G")) # reduced line number from 31 to 21

topHits_drug[topHits_drug$Best_SNP=="rs4410790",]
#--------------------------------------------------------------------------------------------------------
# Get published SNPs for target exposure from GWAS catalog
## URL: https://www.ebi.ac.uk/gwas/search?query=nicotine%20dependence
#--------------------------------------------------------------------------------------------------------
columns_want=c("STRONGEST.SNP.RISK.ALLELE","P.VALUE","OR.or.BETA","MAPPED_TRAIT")

# Get a file with file paths
filePaths=Sys.glob(paste0(locICCRef,"gwas-association-downloaded*.tsv"))

# Reduce input data to one-column SNP lists
for (i in 1:length(filePaths)){
  filePath=filePaths[i]
  # Extract trait name by deleteing date, file prefix and file extension
  ## URL https://stackoverflow.com/questions/43256971/replacing-dates-in-a-character-vector-to-a-specific-format
  yyyy_mm_dd="\\d{2,4}[.-]\\d{2}[.-]\\d{2,4}"
  name1= gsub(basename(filePath),pattern =yyyy_mm_dd,replacement = "" )
  name2=gsub(name1,pattern="gwas-association-downloaded_-",replacement = "")
  name3=gsub(name2,pattern=".tsv",replacement = "")
  trait=gsub(name3,pattern="-",replacement="_")
  
  published_topHits=read.table(filePath,sep="\t",header=T,stringsAsFactors = F) %>% 
                      select_(.dots=columns_want)
  
  ## Split column STRONGEST.SNP.RISK.ALLELE into 2 columns by -
  published_topHits$strongest_SNP <- as.data.frame(do.call("rbind",strsplit(published_topHits$STRONGEST.SNP.RISK.ALLELE,"-"))
                                            ,stringsAsFactors =FALSE )[,1]
  
  published_topHits$risk_allele <- as.data.frame(do.call("rbind",strsplit(published_topHits$STRONGEST.SNP.RISK.ALLELE,"-"))
                                          ,stringsAsFactors =FALSE )[,2] 
  
  # Reduce SNPs to a list of unique SNPs
  published_topHits <- published_topHits %>%
                          filter(P.VALUE < 5e-08) %>%
                            filter(risk_allele !="?") %>% 
                              distinct(strongest_SNP) # 12 obs. of  1 variables
  
  # Name the SNP list by trait
  assign(paste0("published_topSNPs_",trait),published_topHits)
}

#--------------------------------------------------------------
# Join GSCAN QCed GWAS with downloaded GWAS by SNPs
#--------------------------------------------------------------
CYP1A2 <- ImportATabSeparatedFile(input.file.path = paste0(loc.GWAS.catalog,"gwas-association-downloaded_2019-10-21-ensemblMappedGenes_CYP1A2.tsv")
                        ,data.name = "CYP1A2") %>%
  dplyr::select_(.dots=columns_want) %>%
  # Split STRONGEST.SNP.RISK.ALLELE column by dash into 2 columns strongest_SNP, risk_allele 
  tidyr::separate(STRONGEST.SNP.RISK.ALLELE
                  ,c("strongest_SNP","risk_allele")
                  ,sep= "-"
                  ,remove=FALSE) %>%
  dplyr::rename(SNP=strongest_SNP)# dim(CYP1A2) 44 6

# QCed GSCAN GWAS for SI
ImportATabSeparatedFile(input.file.path = paste0(locGSCAN.QC3,"si_noICC.ambiguousSNPRemoved")
                        ,data.name = "GSCAN.GWAS.SI") # dim(GSCAN.GWAS.SI) 12864280       12

GSCAN.GWAS.SI.top.SNPs <- GSCAN.GWAS.SI %>%
  dplyr::filter(PVALUE < 5e-8) %>% 
  dplyr::rename(SNP=RSID)# dim(GSCAN.GWAS.SI.top.SNPs) 238 12

# Clumped GSCAN GWAS for SI
ImportASpaceSeparatedFile(input.file.path = paste0(locGSCAN.QC4,"GWAS_from-clumped-SNPs_si_noICC_LDWindow-kb-10000_R2-0.01_p1-5e-8_p2-1e-6")
                          ,data.name = "clumped.GSCAN.GWAS.SI") # dim(clumped.GSCAN.GWAS.SI) 8 13

# Full join the 3 datasets
full.join <- list( CYP1A2
                  ,GSCAN.GWAS.SI.top.SNPs[,c("SNP","PVALUE")] # .x
                  ,clumped.GSCAN.GWAS.SI[,c("SNP","PVALUE")]) %>% #.y
  purrr::reduce(dplyr::full_join, by= c("SNP")) # dim(full.join) 282 29

rP.GSCAN.PRS.targ.pheno.gp2.wide <- list(rP.GSCAN.PRS.targ.pheno.gp2.wide.allSexes
                                         ,rP.GSCAN.PRS.targ.pheno.gp2.wide.females
                                         ,rP.GSCAN.PRS.targ.pheno.gp2.wide.males) %>% reduce(left_join, by = c("Var1","p.value.threshold"))


