##################################################################################
# Filename: MR_step10-01_function_make-input-data-for-forestplot.R
# Modified from:  /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/forestplot_JueShengOng/JS_forestplot_testscript.R
# Programmer: Chang
# Purpose: (1) Create two output data for forestplot() at step 10-03
# Date created: 20190809
# Reference: https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
# Functions internal: Create.input.data.for.forestplot()
# Note: 
#-----------------------------------------------------------------------------------------

# Type 	File
#------------------------------------------------------------------------------------------------
# Outpu 
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190808  
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

#---------------------------------------------------------------
# Import functions
#---------------------------------------------------------------
#library(Rcpp,lib.loc = "/software/R/R-3.4.1/lib64/R/library")
library(dplyr, lib.loc = "/software/R/R-3.5.1/lib64/R/library")
source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

#----------------------------------------------------------------
# Make a function to create input data for forestplot().
# The function creates two output objects- table data and forest plot data
#----------------------------------------------------------------
#--------------------------------------------
# Testing code with hard-coded arugments
#--------------------------------------------
# input.data.name="outcomes.contin" #"outcomes.binary"
# colname.exposure="exposure"
# colname.exposure.sample="exposure.consortium"
# colname.outcome="outcome"
# colname.outcome.sample="outcome.consortium"
# colname.group1.nsnp="nsnp.IVW"
# colname.group1.survived.multiple.testing="survived.multiple.testing.IVW"
# colname.group1.effect.size="effect.size.estimate.IVW"
# colname.group1.effect.LBound="effect.size.LBound.IVW"
# colname.group1.effect.UBound="effect.size.UBound.IVW"
# colname.group2.survived.multiple.testing="survived.multiple.testing.obs"
# colname.group2.effect.size= "effect.size.estimate.obs"
# colname.group2.effect.LBound="effect.size.LBound.obs"
# colname.group2.effect.UBound="effect.size.UBound.obs"
# output.table.data.name="outcomes.contin.table.data" #"outcomes.binary.table.data"
# output.plot.data.name="outcomes.contin.forestplot.data" #"outcomes.binary.forestplot.data"

Create.input.data.for.forestplot <- function( input.data.name
                                              ,colname.exposure="exposure.label"
                                              ,colname.exposure.sample="exposure.consortium"
                                              ,colname.outcome="outcome.label"
                                              ,colname.outcome.sample="outcome.consortium"
                                              ,colname.group1.nsnp="nsnp.IVW"
                                              ,colname.group1.survived.multiple.testing="survived.multiple.testing.IVW"
                                              ,colname.group1.effect.size="effect.size.estimate.IVW"
                                              ,colname.group1.effect.LBound="effect.size.LBound.IVW"
                                              ,colname.group1.effect.UBound="effect.size.UBound.IVW"
                                              ,colname.group2.survived.multiple.testing="survived.multiple.testing.obs"
                                              ,colname.group2.effect.size= "effect.size.estimate.obs"
                                              ,colname.group2.effect.LBound="effect.size.LBound.obs"
                                              ,colname.group2.effect.UBound="effect.size.UBound.obs"
                                              ,output.table.data.name
                                              ,output.plot.data.name){
  
  # Get input data
  data <- get(input.data.name) # dim(data) 59 28
  
  # Create table header text 
  header.exposure.sample <- "Sample"
  header.exposure.trait <- "Exposure"
  header.outcome.sample <- "Sample"
  header.outcome.trait <- "Outcome"
  header.numb.SNPs <- "N SNPs"
  header.effect.size.MR.IVW <-  "(a)MR IVW estimator"
  header.effect.size.ObsAsso <- "(b)Obser. estimator"
  
  # Specify number of decimal places for rounding values
  dec_size <- 2
  
  # Round estimates (OR) from MR IVW
  effect.sizes.MR.IVW.round <-  as.numeric(format(round(data[,colname.group1.effect.size], dec_size)
                                                  ,nsmall = dec_size)) # length(effect.sizes.MR.IVW.round) 24
  
  effect.sizes.MR.IVW.lower.round <- as.numeric(format(round(data[,colname.group1.effect.LBound], dec_size)
                                                       , nsmall = dec_size))
  
  effect.sizes.MR.IVW.upper.round <- as.numeric(format(round(data[,colname.group1.effect.UBound], dec_size)
                                                       , nsmall = dec_size))
  
  # Round estimates (OR) from observational associations
  effect.sizes.Obs.asso.round <-  as.numeric(format(round(data[,colname.group2.effect.size], dec_size)
                                                    , nsmall = dec_size))
  
  effect.sizes.Obs.asso.lower.round <- as.numeric(format(round(data[,colname.group2.effect.LBound], dec_size)
                                                         , nsmall = dec_size))
  
  effect.sizes.Obs.asso.upper.round <- as.numeric(format(round(data[,colname.group2.effect.UBound], dec_size)
                                                         , nsmall = dec_size)) # length(effect.sizes.Obs.asso.upper.round) 24
  
  # Conditionally concatenate effect sizes, lower and upper limits from MR IVW estimators to a vector. Expect no NA in the data
  ## if p values survived multiple testing, new value= OR [lower, upper] *
  ## if p values did not survive multiple testing, new value= OR [lower, upper]
  plot.data.effect.sizes.95CI.MR.IVW <- ifelse(data[,colname.group1.survived.multiple.testing]==TRUE,
                                               paste0(effect.sizes.MR.IVW.round
                                                      , "["
                                                      , effect.sizes.MR.IVW.lower.round
                                                      , ", "
                                                      , effect.sizes.MR.IVW.upper.round
                                                      , "]"
                                                      ," *"),
                                               paste0(effect.sizes.MR.IVW.round
                                                      , "["
                                                      , effect.sizes.MR.IVW.lower.round
                                                      , ", "
                                                      , effect.sizes.MR.IVW.upper.round
                                                      , "]"))
  
  # Conditionally concatenate effect sizes, lower and upper limits from observational associations to a vector
  ## if vector elements are NA, then set the new value to blank
  ## if vector elements are not NA and p values survived multiple testing, then new value= OR [lower, upper] *
  ## if vector elements are not NA and p values did not survive multiple testing, then new value= OR [lower, upper]
  
  plot.data.effect.sizes.95CI.Obs.asso <- ifelse(is.na(effect.sizes.Obs.asso.round),"",
                                                 ifelse(data[,colname.group2.survived.multiple.testing]==TRUE,
                                                        paste0(effect.sizes.Obs.asso.round
                                                               , "["
                                                               , effect.sizes.Obs.asso.lower.round
                                                               , ", "
                                                               , effect.sizes.Obs.asso.upper.round
                                                               , "]"
                                                               , " *"),
                                                        paste0(effect.sizes.Obs.asso.round
                                                               , "["
                                                               , effect.sizes.Obs.asso.lower.round
                                                               , ", "
                                                               , effect.sizes.Obs.asso.upper.round
                                                               , "]")
                                                 ))
  
  # Fill up first row with some NAs, which will be displayed as blanks in the table
  numb.blank.lines <- 3
  plot.data.effect.sizes.MR.IVW.round <- c(rep(NA,numb.blank.lines),effect.sizes.MR.IVW.round) 
  plot.data.effect.sizes.MR.IVW.lower.round <- c(rep(NA,numb.blank.lines),effect.sizes.MR.IVW.lower.round) 
  plot.data.effect.sizes.MR.IVW.upper.round <- c(rep(NA,numb.blank.lines), effect.sizes.MR.IVW.upper.round) 
  
  plot.data.effect.sizes.Obs.asso.round <- c(rep(NA,numb.blank.lines),effect.sizes.Obs.asso.round) 
  plot.data.effect.sizes.Obs.asso.lower.round<-c(rep(NA,numb.blank.lines),effect.sizes.Obs.asso.lower.round) 
  plot.data.effect.sizes.Obs.asso.upper.round<-c(rep(NA,numb.blank.lines), effect.sizes.Obs.asso.upper.round) 
  
  # Group estimate, lower and upper limits into mean, lower and upper. These will be used to make a forestplot
  ## row.names=c(NA,-27L) seems okay with input data with row number different from 27
  forestplot.data <-structure(list(
    mean= cbind(  plot.data.effect.sizes.MR.IVW.round
                 ,plot.data.effect.sizes.Obs.asso.round)
    ,lower=cbind(  plot.data.effect.sizes.MR.IVW.lower.round
                  ,plot.data.effect.sizes.Obs.asso.lower.round)
    ,upper=cbind(  plot.data.effect.sizes.MR.IVW.upper.round
                  ,plot.data.effect.sizes.Obs.asso.upper.round)
    ,.Names = c("mean", "lower", "upper")
    ,row.names = c(NA, -27L) # -27L: -nrow(data.frame) lines. This should be number of rows in the data.frame
    ,class = "data.frame"))
  
  # Assign user-defined name to forest plot data
  assign(output.plot.data.name,forestplot.data, envir = .GlobalEnv)
  
  # Create table headers. The order of headers should match column order in the table body
  table.headers<-c(  header.exposure.sample, header.exposure.trait
                    ,header.outcome.sample,  header.outcome.trait
                    ,header.numb.SNPs
                    ,header.effect.size.MR.IVW, header.effect.size.ObsAsso) # table header row
  
  # Blank lines between headers and table data
  blank.lines.betw.header.and.body <- rep(NA,length(table.headers)) 
  
  # Create table data as a matrix. Number of table columns should match number of variables in cbind()
  table.body <- cbind(   data[,colname.exposure.sample] 
                        ,data[,colname.exposure] 
                        ,data[,colname.outcome.sample] 
                        ,data[,colname.outcome] 
                        ,data[,colname.group1.nsnp] 
                        ,plot.data.effect.sizes.95CI.MR.IVW
                        ,plot.data.effect.sizes.95CI.Obs.asso) # dim(table.body) 53 7
  
  table.headers.body <- rbind(  table.headers
                               ,blank.lines.betw.header.and.body
                               ,blank.lines.betw.header.and.body
                               ,table.body) # dim(table.headers.body) 56 7
  
  # Assign user-defined object name to the table.headers.body
  assign(output.table.data.name,table.headers.body, envir = .GlobalEnv)
  
} # Close the function

#----------------------------------------------------------------
# An example of calling the function
#----------------------------------------------------------------


#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#