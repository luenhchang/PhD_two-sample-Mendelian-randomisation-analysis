# --------------------------------------------------------------------------------------------------------
# Program           : /mnt/backedup/home/lunC/scripts/MR_ICC_GSCAN_201806/MR_step08-04_heatmap_LDSC-genetic-correlations.R
# Modified from     : 
## /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step22-03_heatmap-genetic-correlations.R 
## /mnt/backedup/home/lunC/scripts/PRS_UKB_201711/PRS_UKB_201711_step18-07-01_heatmap_var-exp-by-PRS_sex-PRS-int-exclu_all-sex-groups.R
# Author            : Chang
# Date created      : 20190814 (Ekka day 2019)
# Purpose           : Plot a heatmap for genetic correlation coefficients from step08-03, showing only significant correlations
# Function external : CalculateCorrBetween2Variables(), CalculateCorrBetween2GroupsOfVariables()
#--------------------------------------------------------------------------------------------------------
# Run dependency    : MR_step08-03_parse-tabulate_LDSC-SNP-heritability_LDSC-genetic-correlations.R

# Type File
#--------------------------------------------------------------------------------------------------------
# Input paste0(loc.LDSC.tabulated,"LDSC-genetic-correlations.tsv")

# Outpu /mnt/backedup/home/lunC/plots/MR_ICC_GSCAN_201806/genetic-correlation-between-use-4-substances.png
#--------------------------------------------------------------------------------------------------------

# Sys.Date()  History
#--------------------------------------------------------------------------------------------------------
# 20191129  Changed multiple testing back to Bonferroni correction. Exported genetic-correlation-between-use-4-substances.png
# 20190928  Changed multiple testing from Bonferroni correction to choose(6,2). Exported genetic-correlation-between-use-4-substances.png 
# 20190924  Reordered row dimension by substance and stages of substance use. Exported genetic-correlation-between-use-4-substances.png
# 20190913  Exported genetic-correlation-between-use-4-substances.png
# 20190814  Exported the 2 PNG files above
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------
# Folder locations under software
#----------------------------------------------------------
loc.R.3.4.1.library <- "/software/R/R-3.4.1/lib64/R/library"
loc.R.3.5.1.library <- "/software/R/R-3.5.1/lib64/R/library"

#----------------------------------------------------------
# Folder locations under my home directory
#----------------------------------------------------------
homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")
locPlots <- paste0(homeDir,"plots/");
loc.plots.manu4 <- paste0(locPlots,"MR_ICC_GSCAN_201806/")
locScripts <- paste0(homeDir,"scripts/MR_ICC_GSCAN_201806/")

#----------------------------------------------------------
# Folder locations under my working directory
#----------------------------------------------------------
workingDir <- "/mnt/lustre/working/lab_nickm/lunC/";
loc.LDSC <- paste0(workingDir,"MR_ICC_GSCAN_201806/LD-score-correlation/")
loc.LDSC.tabulated <- paste0(loc.LDSC,"output/result-tabulated/");

#----------------------------------------------------------
# Import functions
#----------------------------------------------------------
unloadNamespace("tidyselect")
#library(tidyselect,lib.loc = loc.R.3.5.1.library)
library(tidyselect,lib.loc = loc.R.3.4.1.library)
library(dplyr, lib.loc = loc.R.3.4.1.library)
library(corrplot, lib.loc = loc.R.3.4.1.library)
library(RColorBrewer, lib.loc = loc.R.3.4.1.library)

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))
source(paste0(locRFunction,"Rfunction_spread_data_with_duplicated_keys.R"))

#-------------------------------------------------------------------------
# Import tsv data files
## Select rows and columns to create matrixes for use in the heatmap
#-------------------------------------------------------------------------
#columns.want <- c("trait1.consortium.name","trait2.consortium.name","trait1.consortium","trait1.name","trait2.consortium","trait2.name","trait1.consortium.substance","trait2.consortium.substance","rG.esti","rG.p.value","p.value.divided.by.signi.thres")

manu4.rG <- ImportATabSeparatedFile(input.file.path = paste0(loc.LDSC.tabulated,"LDSC-genetic-correlations.tsv")
                                    ,data.name = "manu4.rG") %>% 
  # Create a new column that uniquely identify row dimension of a matrix to create
  dplyr::mutate( trait1.consortium.name= paste(trait1.consortium,trait1.name,sep="_")
                ,trait2.consortium.name= paste(trait2.consortium,trait2.name,sep="_")) %>%
  # Sort data so only half of the off-diaganol parts of the matrix will be populated
  dplyr::arrange(trait1.consortium.name,trait2.consortium.name) %>%
  # Keep data to what are needed
  #dplyr::select_(.dots = columns.want) 
  dplyr::select(trait1.consortium.name,trait2.consortium.name,trait1.consortium,trait1.name,trait2.consortium,trait2.name,trait1.consortium.substance,trait2.consortium.substance,rG.esti,rG.p.value,p.value.divided.by.signi.thres) # dim(manu4.rG) 45 11

# Note the 45 rows will only fill up the lower/upper triangular part of the matrices

#--------------------------------------------------------------------------------------
# Create 3 matrixes with rG, p values, p values divided by significance threshod
#--------------------------------------------------------------------------------------
# Create three empty matrices
## Each matrix should be symmetric and have the same dimension of 10*10
## DON'T use unique levels of either trait1.consortium.name or trait2.consortium.name as the matrix dimensions. You will miss 1. Use both of them
unique(c(manu4.rG$trait1.consortium.name,manu4.rG$trait2.consortium.name))
# [1] "GSCAN_AI"     "GSCAN_CPD"    "GSCAN_DPW"    "GSCAN_SC"     "GSCAN_SI"     "ICC_CI"       "UKB_caffeine"
# [8] "UKB_CPD"      "UKB_ESDPW"    "UKB_PYOS" 

## Order variable by substance (alcohol> caffeine> cannabis> nicotine), and by stages of use (initiation> quantitative use> cessation) 
row.dimension <- c( "GSCAN_DPW","UKB_ESDPW"
                   ,"UKB_caffeine"
                   ,"ICC_CI"
                   ,paste0("GSCAN_",c("SI","AI","CPD")),"UKB_CPD","UKB_PYOS","GSCAN_SC") # length(row.dimension) 10

col.dimension <- row.dimension

# Create 3 matrices with 1s for populating holding rG and p values
## 1s are used as values on the diagnoal, because self-correlation=1 between same variables but this wasn't previously done
mat.rG <- matrix(NA
                 ,nrow=length(row.dimension)
                 ,ncol=length(col.dimension)
                 ,dimnames = list(row.dimension,col.dimension)
                 ) # dim(mat.rG) 10 10
mat.pvalue <- mat.rG

mat.quotient <- mat.rG

# Populate the 3 matrices based on element names from the data.frame
## Ref: https://stackoverflow.com/questions/37931118/r-fill-matrix-based-on-element-names-from-data-frame
for (i in 1:nrow(manu4.rG)){
  trait1 <- manu4.rG[i,"trait1.consortium.name"]
  trait2 <- manu4.rG[i,"trait2.consortium.name"]

  # Populate half of the off-diagonal parts of the genetic correlation matrix by
  ## Matching column 1 value of data.frame and matrix row coordinate
  ## Matching column 2 value of data.frame and matrix column coordinate
  mat.rG[ rownames(mat.rG)== manu4.rG$trait1.consortium.name[i]
         ,colnames(mat.rG)== manu4.rG$trait2.consortium.name[i]] <- manu4.rG[i,"rG.esti"]
  
  # Populate the other half of the off-diagonal parts of the genetic correlation matrix by
  ## Matching column 1 value of data.frame and matrix column coordinate
  ## Matching column 2 value of data.frame and matrix row coordinate
  mat.rG[ rownames(mat.rG)== manu4.rG$trait2.consortium.name[i]
          ,colnames(mat.rG)== manu4.rG$trait1.consortium.name[i]] <- manu4.rG[i,"rG.esti"]
  
  
  # Populate half of the off-diagonal parts of the p value matrix
  mat.pvalue[ rownames(mat.pvalue)== manu4.rG$trait1.consortium.name[i]
             ,colnames(mat.pvalue)== manu4.rG$trait2.consortium.name[i] ] <- manu4.rG[i,"rG.p.value"] 
  
  # Populate the other half of the off-diagonal parts of the p value matrix
  mat.pvalue[ rownames(mat.pvalue)== manu4.rG$trait2.consortium.name[i]
              ,colnames(mat.pvalue)== manu4.rG$trait1.consortium.name[i] ] <- manu4.rG[i,"rG.p.value"] 
  
  # Populate half of the off-diagonal parts of the quotient matrix
  mat.quotient[ rownames(mat.quotient)== manu4.rG$trait1.consortium.name[i]
               ,colnames(mat.quotient)== manu4.rG$trait2.consortium.name[i] ] <- manu4.rG[i,"p.value.divided.by.signi.thres"]
  
  # Populate the other half of the off-diagonal parts of the quotient matrix
  mat.quotient[ rownames(mat.quotient)== manu4.rG$trait2.consortium.name[i]
                ,colnames(mat.quotient)== manu4.rG$trait1.consortium.name[i] ] <- manu4.rG[i,"p.value.divided.by.signi.thres"]
}

# Fill up the diagonal part with 1 
diag(mat.rG) <- 1
diag(mat.pvalue) <- 1
diag(mat.quotient) <- 1

#-----------------------------------------------------------------
# Make a heatmap
#-----------------------------------------------------------------
# Determine lower and upper limits by the range of the off-diagnoal values
mat.rG.diag.as.NA <- mat.rG
diag(mat.rG.diag.as.NA) <- NA
range(mat.rG.diag.as.NA, na.rm = TRUE) # [1] -0.8205  0.9409
rG.lower.limit <- -0.85
rG.upper.limit <- 0.95

# Set the diagonal to 0 so that the legend scale won't be as high as 1
mat.rG.diag.as.0 <- mat.rG
diag(mat.rG.diag.as.0) <- 0

# Calcuate significance threshold cutoff, if not using the quotient matrix
numb.rG.tests <- nrow(manu4.rG) # 45
plot.significance.cutoff <- 0.05/numb.rG.tests

# Significance threshold by trait groups (this is not in use)
# numb.trait.groups <- length(unique(manu4.rG$trait1.consortium.substance))
# numb.pairwise.combinations <- choose(numb.trait.groups,2) # 15
# plot.significance.cutoff <- 0.05/numb.pairwise.combinations

# Group substances together for displaying them in different colors
## Manuscript Color Target phenotypes 
##--------------------------------------------------------
##  4         black  
##--------------------------------------------------------
textLabelColor.vertical.manu4 <- c(rep("black",times=length(row.dimension)))
textLabelColor.horizont.manu4 <- textLabelColor.vertical.manu4

#-----------------------------------------------------------------------------
# Create gradient colors for corrplot
#-----------------------------------------------------------------------------
# Get gradient of 10 colors from color 1 to color 2
## https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
colfunc <- colorRampPalette(c("red","yellow"))
gradient.10.colors.red.to.yellow <- colfunc(10)

# Take a look at the gradient colors
plot(rep(1,10),col=gradient.10.colors.red.to.yellow,pch=19,cex=3)

manu4.y.axis.title <- ""
manu4.x.axis.title <- ""

# Make output figure file path
output.file.path <- paste0(loc.plots.manu4,"genetic-correlation-between-use-4-substances.png")

png(file=output.file.path
    ,width = 1000
    ,height =1000 )

# Make a plot for the lower triangular of the correlation matrix
## Comment this code chunk out when getting the range of color legend
corrplot(mat.rG.diag.as.0 #pheno_corr_coef_mx_diag_as_0
         , col= colorRampPalette(c("red","yellow"))(200)  # return 200 colors between red and yellow
         , is.corr = F
         , method = "square" # # Display the correlation coefficient as the method
         , addCoef.col = "black" # Add correlation coefficients in color black
         , type = "lower"
         , diag = F
         , tl.col="black"
         , cl.cex = 1.25
         , cl.lim = c(rG.lower.limit,rG.upper.limit) # The limits (x1, x2) in the 'c'olor'l'abel.
         #, cl.pos= "n" # position of color legend (labels). cl.pos="n" draws no colorlegend
         , tl.cex = 1.25
         #, tl.pos="n" # rid of labels
         , tl.srt=45 #Text label color and rotation
         , number.cex = 1.5
         , p.mat= mat.pvalue # Matrix of p-value, if NULL, arguments sig.level, insig, pch, pch.col, pch.cex is invalid.
         , sig.level = plot.significance.cutoff
         , insig = "blank")

dev.off() # End the png()

setwd(loc.plots.manu4)

#file.copy(paste0(loc.plots.manu4,"genetic-correlation-between-use-4-substances.png"),paste0(loc.plots.manu4,"FIG1.png"))

# Call the function for making a heatmap 
# source(paste0(locRFunction,"RFunction_correlation-plot_corrplot-modified.R"))
# 
# CallFunctionPlotHeatmapByCorrplot(output.file.path=output.file.path
#                                   ,plot.width=1000
#                                   ,plot.height=1000
#                                   ,plot.data.matrix= mat.rG.diag.as.0
#                                   ,plot.method = "square" # Display the correlation coefficient as the method
#                                   ,plot.type="lower" # lower if only display lower triangular
#                                   ,correlation.or.not=FALSE # Set to FALSE if plotting a general matrix
#                                   ,cex.number= NULL # Cex parameter to send to the call to text when writing corr coeff into
#                                   ,ordering.method="original"
#                                   ,num.rectangles.drew.hierarchical.cluster=NULL
#                                   ,color=gradient.10.colors.red.to.yellow
#                                   ,text.label.horizontal.color= textLabelColor.horizont.manu4 
#                                   ,text.label.vertical.color= textLabelColor.vertical.manu4 
#                                   ,R2.lowerBound=rG.lower.limit
#                                   ,R2.upperBound=rG.upper.limit
#                                   ,plot.label.cex=1.5
#                                   ,plot.data.signi.matrix= mat.quotient
#                                   ,plot.data.signi.threshold= 1 
#                                   ,colorlegend.position.x=c(1,6)
#                                   ,colorlegend.position.y=c(7,7+2)
#                                   ,colorlegend.cex=1.5
#                                   ,YAxisTitle=manu4.y.axis.title
#                                   ,XAxisTitle=manu4.x.axis.title )

