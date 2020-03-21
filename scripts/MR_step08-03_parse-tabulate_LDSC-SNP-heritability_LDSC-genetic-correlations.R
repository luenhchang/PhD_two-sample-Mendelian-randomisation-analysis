#!/usr/bin/Rscript

#---------------------------------------------------------------------------------------------
# Program       : MR_step08-03_parse-tabulate_LDSC-SNP-heritability_LDSC-genetic-correlations.R
# Modified from : PRS_UKB_201711_step22-03_heatmap-genetic-correlations.R
# Date created  : 20190315
# Purpose       : Parse and tabulate data from LDSC log files
# Note: 
# Reference     : LDSC rg log file explanation: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
#----------------------------------------------------------------------------------------
# Run dependency: MR_step00-02_calculate_cases_controls.R (copy sample prevalence of ICC CI, GSCAN SI, and GSCAN SC)
# Function external: ImportATabSeparatedFile(), ExportFileTabSeparated()

# Type  Files
#----------------------------------------------------------------------------------------------
# Input paste0(loc.LDSC.input,"file-info_munged-QCed-GWASs.tsv")
# Input Sys.glob(paste0(loc.LDSC.h2,"SNP-heritability_*-*.log")) (10 files)
# Input Sys.glob(paste0(loc.LDSC.rG,"genetic-correlation_between_*-*_and_*-*.log")) (45 files)

# Outpu paste0(loc.LDSC.tabulated,"LDSC-SNP-heritability_sample-size-prevalence.tsv")
# Outpu paste0(loc.LDSC.tabulated,"LDSC-genetic-correlations.tsv")
#----------------------------------------------------------------------------------------
# Sys.Date()  History
#----------------------------------------------------------------------------------------
# 20191118  Calculated p values from SNP heritability and SE. Exported paste0(loc.LDSC.tabulated,"LDSC-SNP-heritability_sample-size-prevalence.tsv")
# 20191001 Exported paste0(loc.LDSC.tabulated,"LDSC-SNP-heritability_sample-size-prevalence.tsv")   
# 20190928 Exported paste0(loc.LDSC.tabulated,"LDSC-genetic-correlations.tsv")
# 20190928 Changed multiple testing correction from Bonferroni to choose(6,2)=15 (collapse 10 traits to 6 trait groups, each per substance and sample). 
# 20190814 Exported the 2 files above (Exclude cups of coffee per day (CCPD) from SNP heritability and genetic correlations)
# 20190813 Exported the 2 files above
# 20190414 Exported the 2 files above
# 20190412 Exported the 2 files above
#----------------------------------------------------------------------------------------

#----------------------------------------------------------
# Folder locations under my home directory
#----------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locPlots <- paste0(homeDir,"plots/");
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

#----------------------------------------------------------
# Folder locations under my working directory
#----------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
loc.LDSC <- paste0(workingDir,"MR_ICC_GSCAN_201806/LD-score-correlation/")
loc.LDSC.input <- paste0(loc.LDSC,"input/")
loc.LDSC.rG <- paste0(loc.LDSC,"output/genetic-correlations")
loc.LDSC.h2 <- paste0(loc.LDSC,"output/SNP-heritability")
loc.LDSC.tabulated <- paste0(loc.LDSC,"output/result-tabulated/");
#dir.create(loc.LDSC.tabulated)

#----------------------------------------------------------
# Import functions
#----------------------------------------------------------
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#---------------------------------------------------------------------------------------------------
# Import files
#---------------------------------------------------------------------------------------------------
# Add sample size
## Copy N.continuous to sample.size
## Add up N.cases and N.controls for binary traits as sample.size
## Add sample prevalence of ICC CI, GSCAN SC and GSCAN SI. All these prevalences were obtained thru emails with GSCAN and ICC persons, rather than any R scripts
munged.QCed.GWASs <- ImportATabSeparatedFile(input.file.path = paste0(loc.LDSC.input,"file-info_QCed-GWASs.tsv")
                                             ,data.name = "munged.QCed.GWASs") %>% 
  dplyr::mutate(sample.size=dplyr::case_when(trait.type=="continuous" ~ N.continuous
                                     ,trait.type=="binary" ~ N.cases+ N.controls)) %>%
  dplyr::mutate(sample.prevalence=dplyr::case_when( consortium=="GSCAN" & trait=="SC" ~ 0.659
                                                   ,consortium=="GSCAN" & trait=="SI" ~ 0.539
                                                   ,consortium=="ICC" & trait=="CI" ~ 0.267
                                                   ,TRUE ~ as.numeric(NA)))# dim(munged.QCed.GWASs) 11 19 file-info_munged-QCed-GWASs.tsv

#-----------------------------------------------------------------------------------------------------
# Parse LDSC SNP heritability log files
#-----------------------------------------------------------------------------------------------------
pattern.file.name.h2.log <- glob2rx("SNP-heritability_*-*.log") # "^SNP-heritability_.*-.*\\.log$"

h2.log <- data.frame(file.path=list.files(path=loc.LDSC.h2
                                          ,pattern = pattern.file.name.h2.log
                                          ,full.names = TRUE)
                     ,stringsAsFactors = F) # dim(h2.log) 11 1

h2.log <- h2.log %>% 
  # Exclude CCPD
  dplyr::filter(!grepl("CCPD",file.path)) %>%
  dplyr::mutate(file.name=basename(file.path)
                ,consoritum.trait=gsub(file.name, pattern = "SNP-heritability_|.log",replacement = "")) %>%
  tidyr::separate(col=consoritum.trait
                  ,into=c("consortium","trait")
                  ,sep="-"
                  ,remove=TRUE) # dim(h2.log) 10 4
  
h2.log3 <- dplyr::left_join(h2.log
                            ,munged.QCed.GWASs[,c("consortium","trait","substance","trait.type","trait.definition","data.type")]
                            ,by=c("consortium"="consortium"
                                  ,"trait"="trait")) # dim(h2.log3) 10 8

# Create an empty data frame for appending data from each iteration's data frame
append_h2 <- data.frame()

# Parse SNP heritability log file
## Number of iterations: 11
for (i in 1:nrow(h2.log3)){
  file.path <- h2.log3[i,"file.path"] # length(h2.log) 32
  consortium <- h2.log3[i,"consortium"]
  trait <- h2.log3[i,"trait"]
  trait.type <- h2.log3[i,"trait.type"]
  trait.definition <- h2.log3[i,"trait.definition"]
  data.type <- h2.log3[i,"data.type"]
  
  # Read all the lines of a log file
  h2.log.all.lines <- readLines(file.path) # length(h2.log.all.lines) 34
  
  # Find line number with patterns
  pattern.h2.continuous <- "Total Observed scale h2: "
  pattern.h2.binary <- "Total Liability scale h2: "
  
  if (trait.type=="binary"){
    line.numb.h2 <- which(grepl(pattern.h2.binary,h2.log.all.lines)) # 28
  } else if (trait.type=="continuous"){
    line.numb.h2 <- which(grepl(pattern.h2.continuous,h2.log.all.lines)) # 26
  }
  # Read lines by starting line number and length
  h2.log.5lines <- scan(file.path
                        , what=character()
                        , sep='\n'
                        , skip=line.numb.h2-1 # line number to skip before beginning to read data values.
                        , nlines=5 # how many lines to read from the starting line
  )
  
  # Extract values from the lines that are just read
  h2_and_se <- gsub(h2.log.5lines[1],pattern = "Total Observed scale h2: |Total Liability scale h2: |\\(|\\)",replacement="") # "0.0309 0.0054"
  h2 <- as.numeric(unlist(strsplit(h2_and_se,split = " "))[1])
  h2.se <- as.numeric(unlist(strsplit(h2_and_se,split = " "))[2])
  
  lambda.GC <- as.numeric(gsub(h2.log.5lines[2],pattern ="Lambda GC: " ,replacement = ""))
  mean.chisq <- as.numeric(gsub(h2.log.5lines[3],pattern ="Mean Chi\\^2: " ,replacement = ""))
  
  intercept.and.se <- gsub(h2.log.5lines[4],pattern = "Intercept: |\\(|\\)",replacement = "")
  intercept <- as.numeric(unlist(strsplit(intercept.and.se,split = " "))[1])
  intercept.se <- as.numeric(unlist(strsplit(intercept.and.se,split = " "))[2]) # standard error of intercept
  
  ratio <- gsub(h2.log.5lines[5],pattern = "Ratio | \\(usually indicates GC correction\\).",replacement = "")
  
  # Store vectors in a data.frame
  df <- data.frame(consortium=consortium
                   ,trait=trait
                   ,trait.definition=trait.definition
                   ,trait.type=trait.type
                   ,data.type=data.type
                   ,h2=h2
                   ,h2.se=h2.se
                   ,lambda.GC=lambda.GC
                   ,mean.chisq=mean.chisq
                   ,intercept=intercept
                   ,intercept.se=intercept.se
                   ,ratio=ratio
                   ,stringsAsFactors = F)
  
  # Append current iteration to the append data frame
  append_h2 <- rbind(append_h2,df)
}  
# dim(append_h2) 10 12

# Merge sample sizes and SNP heritability
h2.sample.sizes <- dplyr::left_join(append_h2
                                    ,munged.QCed.GWASs[,c("consortium","trait","substance","trait.type","population.prevalence","sample.prevalence","sample.size")]
                                    ,by=c("consortium"="consortium"
                                          ,"trait"="trait"
                                          ,"trait.type"="trait.type")) # dim(h2.sample.sizes) 10 16

h2.sample.sizes2 <- h2.sample.sizes %>% 
  # Calculate Z scores by dividing SNP-h2 by SE
  dplyr::mutate( z=h2/h2.se
  # Calculate two tailed p value from z scores, assuming that SNP-h2/SE follows a Z distribution
                ,pvalue = 2*pnorm(abs(z), lower.tail = F))


# Export data as a TSV
ExportFileTabSeparated(data = h2.sample.sizes2
                       , missing.values.as = "NA"
                       , output.file.path = paste0(loc.LDSC.tabulated,"LDSC-SNP-heritability_sample-size-prevalence.tsv"))

#-----------------------------------------------------------------------------------------------------
# Parse LDSC genetic correlation log files
#-----------------------------------------------------------------------------------------------------
pattern.file.name.rG.log <- glob2rx("genetic-correlation_between_*-*_and_*-*.log") # "^genetic-correlation_between_.*-.*_and_.*-.*\\.log$"

rG.log <- data.frame(file.path=list.files(path=loc.LDSC.rG
                                          ,pattern=pattern.file.name.rG.log
                                          ,full.names = TRUE)
                     ,stringsAsFactors = F) # dim(rG.log) 55 1

rG.log2 <- rG.log %>% 
  # Subset groups. Exclude CCPD from trait 1 and trait 2. This removes 10 from 55 trait pairs
  dplyr::filter(!grepl("CCPD",file.path)) %>% 
  dplyr::mutate(file.name=basename(file.path)) %>%
  tidyr::separate(col=file.name
                  ,into=c(paste0("file.name.part.",c(1:5)))
                  ,sep="_"
                  ,remove=TRUE) %>%
  tidyr::separate(col=file.name.part.3
                  ,into=c(paste0("trait1.",c("consortium","name")))
                  ,sep="-"
                  ,remove=TRUE) %>%
  tidyr::separate(col=file.name.part.5
                  ,into=c(paste0("trait2.",c("consortium","name")))
                  ,sep="-"
                  ,remove=TRUE) %>%
  tidyr::separate(col=trait2.name
                  ,into=c("trait2.name","file.extension")
                  ,sep="\\."
                  ,remove=TRUE) %>%
  dplyr::select_(.dots=c("file.path","trait1.consortium","trait1.name","trait2.consortium","trait2.name"))  # dim(rG.log2) 45 5

# Add grouping variables from munged.QCed.GWASs to genetic correlation data
rG.log3 <- dplyr::left_join(rG.log2
                            ,munged.QCed.GWASs[,c("consortium","trait","substance","trait.type")]
                            ,by=c("trait1.consortium" = "consortium"
                                  ,"trait1.name" = "trait")) # dim(rG.log3) 45 7

rG.log4 <- dplyr::left_join(rG.log3
                            ,munged.QCed.GWASs[,c("consortium","trait","substance","trait.type")]
                            ,by=c("trait2.consortium" = "consortium"
                                  ,"trait2.name" = "trait")) # dim(rG.log4) 45 9

rG.log4 <- data.table::setnames(rG.log4
                     , old=c("substance.x","trait.type.x","substance.y","trait.type.y")
                     , new=c(paste0("trait1.",c("substance","type")),paste0("trait2.",c("substance","type")))) # dim(rG.log4) 45 9

rG.log4 <- rG.log4 %>%
  dplyr::mutate( trait1.consortium.substance=paste0(trait1.consortium,"_",trait1.substance)
                ,trait2.consortium.substance=paste0(trait2.consortium,"_",trait2.substance))

# Find number of genetic correlation tests. This will be used for creating a significance threshold in the loop using Bonferronic correction
#numb.rG.tests <- nrow(rG.log4) # 45
numb.sample.substance.group <- length(unique(rG.log4$trait1.consortium.substance)) # 6
numb.pairwise.groups <- choose(numb.sample.substance.group,2)

# Creat an empty data.frame for appending data from all the iterations
append.rG <- data.frame()

for (i in 1:nrow(rG.log4)){
  # Extract elements from the file to loop thru
  file.path <- rG.log4[i,"file.path"]
  trait1.consortium <- rG.log4[i,"trait1.consortium"]
  trait1.name <- rG.log4[i,"trait1.name"]
  trait1.substance <- rG.log4[i,"trait1.substance"]
  trait1.type <- rG.log4[i,"trait1.type"]
  trait1.consortium.substance <- rG.log4[i,"trait1.consortium.substance"]
  
  trait2.consortium <- rG.log4[i,"trait2.consortium"]
  trait2.name <- rG.log4[i,"trait2.name"]
  trait2.substance <- rG.log4[i,"trait2.substance"]
  trait2.type <- rG.log4[i,"trait2.type"]
  trait2.consortium.substance <- rG.log4[i,"trait2.consortium.substance"]
  
  # Read all lines of the file
  rG.log.all.lines <- readLines(file.path) # length(rG.log.all.lines) 67
  
  # Extract line number with Genetic correlation summary
  line.numb.rG <- which(grepl("Genetic Correlation: ",rG.log.all.lines)) #57
  
  # Read selective lines of the file
  rG.log.3lines <- scan(file.path
                        , what=character()
                        , sep='\n'
                        , skip=line.numb.rG-1 
                        , nlines=3)
  
  rG.and.SE <- gsub(rG.log.3lines[1],pattern="Genetic Correlation: |\\(|\\)",replacement = "")
  rG.esti <- as.numeric(strsplit(rG.and.SE," ")[[1]][1])
  rG.SE <- as.numeric(strsplit(rG.and.SE," ")[[1]][2])
  rG.Z.score <- as.numeric(gsub(rG.log.3lines[2],pattern="Z-score: ",replacement = ""))
  rG.p.value <- as.numeric(gsub(rG.log.3lines[3],pattern="P: ",replacement = ""))
  
  # Add extracted values to a data.frame
  df.here <- data.frame(trait1.consortium=trait1.consortium
                        ,trait1.name=trait1.name
                        ,trait1.substance=trait1.substance
                        ,trait1.type=trait1.type
                        ,trait1.consortium.substance=trait1.consortium.substance
                        ,trait2.consortium=trait2.consortium
                        ,trait2.name=trait2.name
                        ,trait2.substance=trait2.substance
                        ,trait2.type=trait2.type
                        ,trait2.consortium.substance=trait2.consortium.substance
                        ,rG.esti=rG.esti
                        ,rG.SE=rG.SE
                        ,rG.Z.score=rG.Z.score
                        ,rG.p.value=rG.p.value
                        ,stringsAsFactors = F) %>%
    # Add Bonferroni correction significance threshold 
    # Calculate quotient of p value divided by the threshold. value < 1 means significant, >= 1 means not
    # Indicate whether a LDSC rG test survived multiple testing
    dplyr::mutate(signi.threshold= 0.05/numb.pairwise.groups
                  ,p.value.divided.by.signi.thres= rG.p.value/signi.threshold
                  ,significant_or_not=dplyr::case_when(rG.p.value < signi.threshold ~ "yes"
                                         , TRUE ~ "no")) # dim(df.here) 1 14
  # Append current iteration result to the append data frame
  append.rG <- rbind(df.here,append.rG)
}

dim(append.rG) # 45 17

# Count number of trait pairs survived multiple testing
append.rG.significant <- append.rG %>% dplyr::filter(significant_or_not=="yes") # nrow(append.rG.significant) 33

# Export data as a TSV
ExportFileTabSeparated(data = append.rG
                       , missing.values.as = "NA"
                       , output.file.path = paste0(loc.LDSC.tabulated,"LDSC-genetic-correlations.tsv"))

setwd(locScripts)
#file.copy("MR_step08-03_parse-tabulate_LDSC-SNP-heritability_LDSC-genetic-correlations.R","MR_step09-01_make-table1-data-sources-sample-sizes-types-of-analysis-performed-on-substance-use-traits.R") 
