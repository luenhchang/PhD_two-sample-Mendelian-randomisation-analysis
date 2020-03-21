##################################################################################
# Filename: MR_step10-02_function_create-table-forestplot.R
# Modified from:  /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/forestplot_JueShengOng/JS_forestplot_testscript.R
# Programmer: Chang
# Purpose: Create a function to make a forest plot using forestplot()
# Date created: 20190809
# Dependent:  /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step10-01_function_make-input-data-for-forestplot.R
# Reference: https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
# Functions internal: Make.table.forestplot()
# Note: 
#-----------------------------------------------------------------------------------------

# Type 	File
#------------------------------------------------------------------------------------------------
# Input	paste0(loc.obs.assoc,"results_binary-logistic-regression_linear-regression_full-parameters.tsv")
# Input	paste0(loc.twoSampleMR.tabulated,"MR-analysis-results_all-trait-pairs.tsv")

# Outpu 
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190806  
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Folder locations under my home
#-----------------------------------------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
loc.plot <- paste0(homeDir,"plots/MR_ICC_GSCAN_201806/")
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")
loc.test <- paste0(homeDir,"test/")
dir_ukbPheno <- paste0(homeDir,"data/UKBiobank_phenotype/")

#-----------------------------------------------------------------------------------------
# Folder locations under my working directory
#-----------------------------------------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
locMR <- paste0(workingDir,"MR_ICC_GSCAN_201806/") # location of outcome data
loc.obs.assoc <- paste0(locMR,"observational-associations/")
loc.twoSampleMR.tabulated <- paste0(locMR,"two-sample-MR/result-tabulated")

#----------------------------------------------------------------------------
# Import functions
#----------------------------------------------------------------------------
#library(Rcpp,lib.loc = "/software/R/R-3.4.1/lib64/R/library")
#library(dplyr, lib.loc = "/software/R/R-3.4.1/lib64/R/library")
library(forestplot, lib.loc ="/software/R/R-3.4.1/lib64/R/library")
source(paste0(locScripts,"MR_step10-01_function_make-input-data-for-forestplot.R"))

#----------------------------------------------------------------------------
# Make a function to create a single table and forestplot in 1-by-1 dimension
#----------------------------------------------------------------------------
# Testing code with hard-coded function arguments
# forestplot.width=2500
# forestplot.height=1000
# output.file.path=output.file.path.binary.outcomes
# plot.title=plot.title.outcome.binary
# x.axis.title.text="Odds Ratio [95% CI]"
# input.table.data.name="outcomes.binary.table.data"
# input.forestplot.data.name="outcomes.binary.forestplot.data"
# width.hori.line.columns=7
# table.font.size.cex=1.5
# x.axis.upper.limit=3
# x.axis.lower.limit=-20
# symbol.size=0.35

Make.table.forestplot <- function(forestplot.width=1100
                                  ,forestplot.height=600
                                  ,output.file.path=output.file.path.binary.outcomes
                                  ,plot.title=plot.title.outcome.binary
                                  ,x.axis.title.text="Odds Ratio [95% CI]"
                                  ,input.table.data.name="outcomes.binary.table.headers.body"
                                  ,input.forestplot.data.name="plot.data.outc.binary"
                                  ,width.hori.line.columns=7
                                  ,table.font.size.cex=1.5
                                  ,x.axis.upper.limit=3
                                  ,x.axis.lower.limit=-20
                                  ,symbol.size=0.35){
  
  # Get input data
  table.data <- get(input.table.data.name) # dim(table.data) 34 7
  plot.data <- get(input.forestplot.data.name) # length(plot.data) a list of 6  
  
  # Determine the range of ticks on x axis
  plot.value.ranges <- range(as.vector(plot.data$mean)
                             ,as.vector(plot.data$lower)
                             ,as.vector(plot.data$upper)
                             ,na.rm = TRUE)
  # Set lower limit for x axis
  if (plot.value.ranges[1] <= x.axis.lower.limit ){
    x.axis.tick.minimum <- x.axis.lower.limit
  }
  else {
    x.axis.tick.minimum <- plot.value.ranges[1] - 0.1
  }
  
  
  # Set upper limit for x axis
  if (plot.value.ranges[2] <= x.axis.upper.limit ){
    x.axis.tick.maximum <- plot.value.ranges[2] + 0.1
  }
  else {
    x.axis.tick.maximum <- x.axis.upper.limit
  }
  
  # Gives 10 equally spaced numbers from 0 to 1 
  x.axis.ticks.round <- round(seq( from=x.axis.tick.minimum
                                   ,to=x.axis.tick.maximum
                                   ,length.out = 10)
                              ,1)
  
  numb.col.plot.dim.x <- 1
  numb.col.plot.dim.y <- 1
  
  plot.width <- numb.col.plot.dim.x*forestplot.width
  plot.height <- numb.col.plot.dim.y*forestplot.height
  
  png(file= output.file.path
      ,width = plot.width
      ,height=  plot.height )
  
  par(mfrow=c(numb.col.plot.dim.y,numb.col.plot.dim.x)
      ,mar=c(8, 15, 1, 0.5) # margin of subplots, where x and y axis titles are shown
      ,oma=c(1.5, 2, 1, 1) # outer margin of entire plot
  )
  
  plot1 <- forestplot(labeltext=table.data
                      , mar = unit(rep(20, times = 4), "mm")
                      ,legend = c("(a) MR estimate", "(b) Observational estimate")
                      ,hrzl_lines = list("3"=gpar(lwd=2
                                                  , columns=1:length(width.hori.line.columns))) # width of the horizontal line between table header and body
                      ,colgap=unit(3,"mm") # gap between columns, defaults to 6 mm but for relative widths
                      ,mean= plot.data$mean 
                      ,lower=plot.data$lower 
                      ,upper=plot.data$upper 
                      ,boxsize= symbol.size # size of estimate symbols
                      ,lwd.ci=1.75 # line width (lwd) for the confidence bands
                      ,xticks= x.axis.ticks.round 
                      ,clip=c(0.6,1.4)
                      ,title= plot.title
                      ,zero = 1
                      ,vertices=TRUE
                      ,lineheight = "auto"
                      ,xlab = x.axis.title.text
                      # Change font settings in the table region
                      ,txt_gp = fpTxtGp(label = gpar(fontfamily = "Arial")
                                        ,ticks = gpar(cex=table.font.size.cex) # 1.25
                                        ,xlab  = gpar(cex=table.font.size.cex) # X axis title size # 1.5
                                        ,cex= table.font.size.cex # The font size
                      )
                      ,col=fpColors(box=c("royalblue", "red"),
                                    line=c("darkblue", "red"),
                                    summary=c("darkblue", "red")),
                      new_page = FALSE)
  dev.off()
}

#----------------------------------------------------------------
# An example of calling the function
#----------------------------------------------------------------
# # Make a plot title
# plot.title.outcome.binary <- "Comparison of observational and MR estimates for the association between exposure and binary outcomes"
# 
# Create.input.data.for.forestplot(input.data.name="outcomes.binary"
#                                  ,output.table.data.name="outcomes.binary.table.data"
#                                  ,output.plot.data.name="outcomes.binary.forestplot.data")
# # Set output file path
# output.file.path.binary.outcomes <- paste0(loc.test, "manu4_odds-ratio-95CI_observational-association_MR-IVW.png")
# 

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#