##################################################################################
# Filename: MR_step10-03_forest-plot_effect-size-estimates-95percent-CI_observational-associations_MR-IVW.R
# Modified from:  /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/forestplot_JueShengOng/JS_forestplot_testscript.R
# Programmer: Chang
# Purpose: (1) (2)
# Date created: 20190806
# Reference: https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
# Dependent: /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step10-02_function_create-table-forestplot.R
# Note: 
#-----------------------------------------------------------------------------------------

# Type 	File
#------------------------------------------------------------------------------------------------
# Input	paste0(loc.obs.assoc,"results_binary-logistic-regression_linear-regression_full-parameters.tsv")
# Input	paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")

# Outpu paste0(loc.plot, "manu4_converted-effect-size-95CI_MR-IVW_observational-association_outcomes-*.png")

#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190910  Exported 5 plots above (one plot per outcome type and outcome substance)
# 20190821  Exported the 2 png files above  
# 20190814  Exported the 2 png files above  
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Folder locations under my home
#-----------------------------------------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
loc.plot <- paste0(homeDir,"plots/MR_ICC_GSCAN_201806/")
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
dir_ukbPheno <- paste0(homeDir,"data/UKBiobank_phenotype/")

#-----------------------------------------------------------------------------------------
# Folder locations under my working directory
#-----------------------------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/") # location of outcome data
loc.obs.assoc <- paste0(locMR,"observational-associations/")
loc.twoSampleMR.tabulated <- paste0(locMR,"two-sample-MR/result-tabulated")

#-----------------------------------------------------------------------------------------
# Import functions
#-----------------------------------------------------------------------------------------
#library(Rcpp,lib.loc = "/software/R/R-3.5.1/lib64/R/library")
# Detach all loaded packages
library(dplyr, lib.loc = "/software/R/R-3.4.1/lib64/R/library")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))
source(paste0(locScripts,"MR_step10-02_function_create-table-forestplot.R"))

#-----------------------------------------------------------------------------------------
# Manipulate input data
#-----------------------------------------------------------------------------------------
# Specify common column names in desired order from multiple imported data. Note names are copied from the step belows
columns.common <- c("type.analysis"
                    ,"var.type.outcome","outcome.consortium","outcome.substance","outcome"
                    ,"var.type.exposure","exposure.consortium","exposure.substance","exposure"
                    ,"nsnp"
                    ,"p.value","signi.threshold","survived.multiple.testing"
                    ,"effect.size.estimate","effect.size.LBound","effect.size.UBound") # length(columns.common)
 
#-----------------------------------------------------------------------#
# Import observation association results from selected predictors of 4 substances- coffee, alcohol, smoking and cannabis in the UK Biobank
#-----------------------------------------------------------------------#
significance.threshold.obs.asso <- 0.05/5 # only use first 5 iterations
  
manu4.obser.assoc.m1to5 <- ImportATabSeparatedFile(input.file.path = paste0(loc.obs.assoc,"results_binary-logistic-regression_linear-regression_full-parameters.tsv")
                                                   ,data.name = "manu4.observational.assoc") %>%
  dplyr::filter(iteration %in% c(1:5) & predictor.label %in% c("ECCPD","ESDPW","PYOS","CI","SI")) %>% 
  # Exclude predictor=pack years of smoking and outcome= ever smoked (SI). PYOS is derived from SI.
  dplyr::filter(!(predictor.label == "PYOS" & dep.var.label == "SI")) %>%
  
  # Group outcome and predictor traits to 4 substances
  dplyr::mutate(exposure.substance=case_when( predictor.label== "ECCPD" ~ "caffeine"
                                             ,predictor.label== "ESDPW" ~ "alcohol"
                                             ,predictor.label %in% c("PYOS","SI") ~ "nicotine"
                                             ,predictor.label== "CI" ~ "cannabis" )
                ,outcome.substance=case_when(  dep.var.label== "ECCPD" ~ "caffeine"
                                             ,dep.var.label== "ESDPW" ~ "alcohol"
                                             ,dep.var.label %in% c("PYOS","SI") ~ "nicotine"
                                             ,dep.var.label== "CI" ~ "cannabis" )
                # Add additional grouping variables
                ,exposure.consortium="UKB"
                ,outcome.consortium="UKB"
                ,type.analysis="observational.association"
                ,nsnp=NA
                ,signi.threshold=significance.threshold.obs.asso) %>%

  # Check whether p.value is significant or non-significant
  dplyr::mutate(survived.multiple.testing=case_when( p.value < signi.threshold ~ TRUE
                                                    ,p.value >= signi.threshold ~ FALSE)) %>%
  
  # Rename columns to same colnames as MR.IVW data for merging
  dplyr::rename(exposure = predictor.label
               ,var.type.exposure= predictor.type
               ,outcome= dep.var.label
               ,var.type.outcome=dep.var.type) %>%
  # Select common columns
  dplyr::select_(.dots = columns.common) # dim(manu4.obser.assoc.m1to5) 19 16 

#-------------------------------------------------
# Import MR IVW results 
#-------------------------------------------------
# All binary trait labels
MR.traits.binary <- c("CI","SC","SI")

# All continuous trait labels
MR.traits.contin <- c("AI","caffeine","CPD","DPW","ESDPW","PYOS")

# Determine significance threshold using Bonferroni correction
numb.traits.MR <- 10 # CPD is in GSCAN and UKB, counted as 2
numb.tests.MR <- choose(numb.traits.MR,2)
numb.clumping.p.values <- length(c(5e-08,1e-05))
numb.MR.directions <- length(c("exposure -> outcome","outcome -> exposure"))

significance.threshold.MR.IVW <- 0.05/(numb.tests.MR*numb.clumping.p.values*numb.MR.directions)

# Import and manipulate data # dim(manu4.MR.IVW) 
manu4.MR.IVW <- ImportATabSeparatedFile(input.file.path = paste0(loc.twoSampleMR.tabulated,"/MR-analysis-results_all-trait-pairs.tsv") 
                                        ,data.name = "manu4.MR.results") %>% 
  # Subset groups
  dplyr::filter(exposure %in% c(MR.traits.binary,MR.traits.contin) & outcome %in% c(MR.traits.binary,MR.traits.contin) & exposure.clumping.p1==5e-08 & method.abbreviated=="IVW") %>% 
  ## Exclude trait pairs with similar exposure and outcome definition
  dplyr::filter(!((exposure=="CPD" & outcome=="CPD")|(exposure=="ESDPW" & outcome=="DPW")|(exposure=="DPW" & outcome=="ESDPW"))) %>%
  dplyr::mutate(type.analysis="MR.IVW") %>%
  
  # Add additional grouping variables
  # Check whether p.value is significant or non-significant
  dplyr::mutate( signi.threshold=significance.threshold.MR.IVW
                 ,survived.multiple.testing=case_when( pval < signi.threshold ~ TRUE
                                                       ,pval >= signi.threshold ~ FALSE)) %>%
  # Rename variables for merging purposes
  dplyr::rename(p.value=pval, SE=se) %>%
  
  # Recode caffeine as ECCPD in exposure and outcome for merging purposes
  dplyr::mutate( exposure= dplyr::recode(exposure, `caffeine`="ECCPD")
                 ,outcome= dplyr::recode(outcome, `caffeine`="ECCPD")) %>%
  dplyr::select_(.dots =  columns.common) # dim(manu4.MR.IVW) 77 16

#-------------------------------------------------------------------------------------#
# Horizontally combine the two data above
#-------------------------------------------------------------------------------------#
# All columns from left table are suffixed with .x and from right table are suffixed with .y
## The resulting data have row number similar to the larger right table. Good.
t <- dplyr::full_join( manu4.obser.assoc.m1to5
                      ,manu4.MR.IVW
                      ,by=c(  "var.type.outcome"="var.type.outcome"
                             ,"outcome"="outcome"
                             ,"outcome.consortium"="outcome.consortium"
                             ,"outcome.substance"="outcome.substance"
                             ,"var.type.exposure"="var.type.exposure"
                             ,"exposure"="exposure"
                             ,"exposure.consortium"="exposure.consortium"
                             ,"exposure.substance"="exposure.substance")) # dim(t) 90 24 # "type.analysis"="type.analysis"

# Replace .x with .obs, .y with .IVW
new.header <- gsub(colnames(t),pattern = "\\.x", replacement="\\.obs") %>% 
  gsub(.,pattern = "\\.y", replacement="\\.IVW")

colnames(t) <- new.header

# Add labels to exposure/predictor and outcome using external object manu4.var.name.label
t1 <- dplyr::left_join(t
                       ,manu4.var.name.label
                       ,by=c("exposure"="var.abbre")) %>%
  dplyr::rename(exposure.label=var.label) %>%
  dplyr::left_join(.
                   ,manu4.var.name.label
                   ,by=c("outcome"="var.abbre")) %>%
  dplyr::rename(outcome.label=var.label) # dim(t1) 90 26

# Split the joined file into binary outcomes and continuous outcomes. They will be plotted in two forestplots. Not in one because of different scales between odds ratio and beta
t1.b <- t1 %>% dplyr::filter(var.type.outcome=="binary")  # dim(outcomes.binary) 31 28
outcomes.binary <- t1.b[order( t1.b$outcome.substance
                          ,t1.b$outcome
                          ,t1.b$exposure.substance
                          ,t1.b$exposure),]

# Split data by outcome substance
outcomes.binary.list <- split(outcomes.binary,outcomes.binary$outcome.substance)

t1.c <- t1 %>% dplyr::filter(var.type.outcome=="continuous") # dim(outcomes.contin) 63 26
outcomes.contin <- t1.c[order( t1.c$outcome.substance
                               ,t1.c$outcome
                               ,t1.c$exposure.substance
                               ,t1.c$exposure),]

# Split data by outcome substance
outcomes.contin.list <- split(outcomes.contin, outcomes.contin$outcome.substance)

#----------------------------------------------------------------
# Make forest plots per binary outcomes cannabis
#----------------------------------------------------------------
# Get data from the list
outcome.binary.cannabis <- outcomes.binary.list[["cannabis"]] # dim(outcome.binary.cannabis) 12 26

# Make a plot title
plot.title.outcome.binary.cannabis <- "Comparison of MR IVW and observational estimates for the association between  any exposure and cannabis initiation"

source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))

Create.input.data.for.forestplot(input.data.name="outcome.binary.cannabis"
                                 ,output.table.data.name="outcomes.binary.cannabis.table.data"
                                 ,output.plot.data.name="outcomes.binary.cannabis.forestplot.data")
# Set output file path
file.name.prefix <- "manu4_converted-effect-size-95CI_MR-IVW_observational-association"
file.name.suffix <- "outcomes-binary-cannabis"

output.file.path.binary.outcomes.cannabis <- paste0(loc.plot
                                                    ,file.name.prefix
                                                    ,"_"
                                                    ,file.name.suffix
                                                    ,".png")
# Make a plot
source(paste0(locScripts,"MR_step10-02_function_create-table-forestplot.R"))

Make.table.forestplot( forestplot.width=2000
                      ,forestplot.height=1500
                      ,output.file.path=output.file.path.binary.outcomes.cannabis
                      ,plot.title=plot.title.outcome.binary.cannabis
                      ,x.axis.title.text="Odds Ratio [95% CI]"
                      ,input.table.data.name="outcomes.binary.cannabis.table.data"
                      ,input.forestplot.data.name="outcomes.binary.cannabis.forestplot.data"
                      ,width.hori.line.columns=6
                      ,table.font.size.cex=2
                      ,x.axis.upper.limit=3
                      ,symbol.size=0.25)

setwd(loc.plot)

#----------------------------------------------------------------
# Make forest plots per binary outcomes nicotine
#----------------------------------------------------------------
# Plot binary outcome nicotine
outcome.binary.nicotine <- outcomes.binary.list[["nicotine"]] # dim(outcome.binary.nicotine) 19 26

# Make a plot title
plot.title.outcome.binary.nicotine <- "Comparison of MR IVW and observational estimates for the association between exposure and binary outcomes related to nicotine use"

source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))

Create.input.data.for.forestplot(input.data.name="outcome.binary.nicotine"
                                 ,output.table.data.name="outcomes.binary.nicotine.table.data"
                                 ,output.plot.data.name="outcomes.binary.nicotine.forestplot.data")
# Set output file path
file.name.prefix <- "manu4_converted-effect-size-95CI_MR-IVW_observational-association"
file.name.suffix <- "outcomes-binary-nicotine"

output.file.path.binary.outcomes.nicotine <- paste0(loc.plot
                                                    ,file.name.prefix
                                                    ,"_"
                                                    ,file.name.suffix
                                                    ,".png")
# Make a plot
source(paste0(locScripts,"MR_step10-02_function_create-table-forestplot.R"))

Make.table.forestplot( forestplot.width=2000
                       ,forestplot.height=1500
                       ,output.file.path=output.file.path.binary.outcomes.nicotine
                       ,plot.title= plot.title.outcome.binary.nicotine
                       ,x.axis.title.text="Odds Ratio [95% CI]"
                       ,input.table.data.name="outcomes.binary.nicotine.table.data"
                       ,input.forestplot.data.name="outcomes.binary.nicotine.forestplot.data"
                       ,width.hori.line.columns=7
                       ,table.font.size.cex=2
                       ,x.axis.upper.limit=3
                       ,symbol.size=0.25)

#----------------------------------------------------------------
# Make a plot for continuous outcomes alcohol
#----------------------------------------------------------------
# Get data from the list
outcomes.contin.alcohol <- outcomes.contin.list[["alcohol"]] # dim(outcomes.contin.alcohol) 16 26

# Make a plot title
plot.title.outcome.contin.alcohol <- "Comparison of MR IVW and observational estimates for the association between exposure and continuous outcomes related to alcohol use"

source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))

Create.input.data.for.forestplot(input.data.name="outcomes.contin.alcohol"
                                 ,output.table.data.name="outcomes.contin.alcohol.table.data"
                                 ,output.plot.data.name="outcomes.contin.alcohol.forestplot.data")
# Set output file path
file.name.prefix <- "manu4_converted-effect-size-95CI_MR-IVW_observational-association"
file.name.suffix <- "outcomes-continous-alcohol"
output.file.path.contin.outcomes.alcohol <- paste0(loc.plot
                                                    ,file.name.prefix
                                                    ,"_"
                                                    ,file.name.suffix
                                                    ,".png")
# Make a plot
Make.table.forestplot( forestplot.width=2500
                      ,forestplot.height=2000
                      ,output.file.path=output.file.path.contin.outcomes.alcohol
                      ,plot.title=plot.title.outcome.contin.alcohol
                      ,x.axis.title.text="Effect size [95% CI]"
                      ,input.table.data.name="outcomes.contin.alcohol.table.data"
                      ,input.forestplot.data.name="outcomes.contin.alcohol.forestplot.data"
                      ,width.hori.line.columns=7
                      ,table.font.size.cex=2.5
                      ,x.axis.upper.limit=3
                      ,symbol.size=0.25)

#----------------------------------------------------------------
# Make a plot for continuous outcomes caffeine
#----------------------------------------------------------------
# Get data from the list
outcomes.contin.caffeine <- outcomes.contin.list[["caffeine"]] # dim(outcomes.contin.caffeine) 10 26

# Make a plot title
plot.title.outcome.contin.caffeine <- "Comparison of MR IVW and observational estimates for the association between exposure and caffeine consumption (outcomes)"

source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))

Create.input.data.for.forestplot(input.data.name="outcomes.contin.caffeine"
                                 ,output.table.data.name="outcomes.contin.caffeine.table.data"
                                 ,output.plot.data.name="outcomes.contin.caffeine.forestplot.data")
# Set output file path
file.name.prefix <- "manu4_converted-effect-size-95CI_MR-IVW_observational-association"
file.name.suffix <- "outcomes-continous-caffeine"
output.file.path.contin.outcomes.caffeine <- paste0(loc.plot
                                                   ,file.name.prefix
                                                   ,"_"
                                                   ,file.name.suffix
                                                   ,".png")
# Make a plot
source(paste0(locScripts,"MR_step10-02_function_create-table-forestplot.R"))

Make.table.forestplot( forestplot.width=2500
                       ,forestplot.height=2000
                       ,output.file.path=output.file.path.contin.outcomes.caffeine
                       ,plot.title=plot.title.outcome.contin.caffeine
                       ,x.axis.title.text="Effect size [95% CI]"
                       ,input.table.data.name="outcomes.contin.caffeine.table.data"
                       ,input.forestplot.data.name="outcomes.contin.caffeine.forestplot.data"
                       ,width.hori.line.columns=7
                       ,table.font.size.cex=3
                       ,x.axis.upper.limit=60
                       ,x.axis.lower.limit=-30
                       ,symbol.size=0.25)

#----------------------------------------------------------------
# Make a plot for continuous outcomes nicotine
#----------------------------------------------------------------
# Get data from the list
outcomes.contin.nicotine <- outcomes.contin.list[["nicotine"]] # dim(outcomes.contin.nicotine) 33 26

# Make a plot title
plot.title.outcome.contin.nicotine <- "Comparison of MR IVW and observational estimates for the association between any exposure and nicotine use (continuous outcomes)"

source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))

Create.input.data.for.forestplot(input.data.name="outcomes.contin.nicotine"
                                 ,output.table.data.name="outcomes.contin.nicotine.table.data"
                                 ,output.plot.data.name="outcomes.contin.nicotine.forestplot.data")
# Set output file path
file.name.prefix <- "manu4_converted-effect-size-95CI_MR-IVW_observational-association"
file.name.suffix <- "outcomes-continous-nicotine"
output.file.path.contin.outcomes.nicotine <- paste0(loc.plot
                                                    ,file.name.prefix
                                                    ,"_"
                                                    ,file.name.suffix
                                                    ,".png")
# Make a plot
Make.table.forestplot( forestplot.width=3500
                       ,forestplot.height=2500
                       ,output.file.path=output.file.path.contin.outcomes.nicotine
                       ,plot.title=plot.title.outcome.contin.nicotine
                       ,x.axis.title.text="Effect size [95% CI]"
                       ,input.table.data.name="outcomes.contin.nicotine.table.data"
                       ,input.forestplot.data.name="outcomes.contin.nicotine.forestplot.data"
                       ,width.hori.line.columns=7
                       ,table.font.size.cex=3
                       ,x.axis.upper.limit=10
                       ,symbol.size=0.5)

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#