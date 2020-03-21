##################################################################################
# filename: MR_step06-04_plot_MR-analysis-results.R
# program author: Chang
# purpose: Plot SNP effect on exposure against SNP effect on outcome using Mendelian randomisation results from MR_step06-03
# date created: 20190325
# Function internal: mr_scatter_plot_modified()
# Function external: ImportATabSeparatedFile()
# file directory: 
#-----------------------------------------------------------------------------------------
# Type 	File
#------------------------------------------------------------------------------------------------
# Input paste0(loc.harmonised.data,"harmonised-data_exposure-UKB-LDClumping-CCPD_outcome-ICC-CI.tsv")
# Input	paste0(loc.harmonised.data,"harmonised-data_exposure-UKB-GCTACOJO-CCPD_outcome-ICC-CI.tsv")

# Outpu paste0(locPlots_MR,"zfig001-01_MR-SNP-effect_exposure-UKB-CCPD-LDClumped_outcome-ICC-cannabis-initiation.png"
# Outpu paste0(locPlots_MR,"zfig001-02_MR-SNP-effect_exposure-UKB-CCPD-GCTACOJO_outcome-ICC-cannabis-initiation.png")
#------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# Sys.time()  Update
#-----------------------------------------------------------------------------------------
# 20190326  Exported the 2 PNG files above
#-----------------------------------------------------------------------------------------

homeDir <- "/mnt/backedup/home/lunC/";
locRFunction <- paste0(homeDir,"scripts/RFunctions/")

workingDir <- "/mnt/lustre/working/lab_nickm/lunC/"
loc.MR <- paste0(workingDir,"MR_ICC_GSCAN_201806/")
loc.GWAS.catalog <- paste0(loc.MR,"GWAS-catalog-download/")
loc.harmonised.data <- paste0(loc.MR,"twoSampleMR-harmonised-data/")

loc.MR.report <- paste0(workingDir,"MR_ICC_GSCAN_201806/result_MR_reports/")
loc.MR.report.GSCAN.LDClumping <- paste0(loc.MR.report,"GSCAN-LDClumping")
loc.MR.report.UKB.LDClumping <- paste0(loc.MR.report,"UKB-LDClumping")
loc.MR.report.UKB.GCTACOJO <- paste0(loc.MR.report,"UKB-GCTACOJO")

source(paste0(locRFunction,"RFunction_import_export_single_file.R"))
source(paste0(locRFunction,"RFunction_format-values.R"))

# Import files
ImportATabSeparatedFile(input.file.path = paste0(loc.harmonised.data,"harmonised-data_exposure-UKB-LDClumping-CCPD_outcome-ICC-CI.tsv")
                        ,data.name = "harmonised.CI.UKB.LDClumped.CCPD")

ImportATabSeparatedFile(input.file.path = paste0(loc.harmonised.data,"harmonised-data_exposure-UKB-GCTACOJO-CCPD_outcome-ICC-CI.tsv")
                        ,data.name = "harmonised.CI.UKB.GCTACOJO.CCPD")

ImportATabSeparatedFile(input.file.path = paste0(loc.GWAS.catalog,"EFO-0004330_coffee-consumption_downloaded-20190327.tsv")
                        ,data.name = "GWAS.catalog.EFO.0004330") # dim(GWAS.catalog.EFO.0004330) 34 38

ImportATabSeparatedFile(input.file.path = paste0(loc.GWAS.catalog,"EFO-0007585_cannabis-use_downloaded-20190328.tsv")
                        ,data.name = "GWAS.catalog.EFO.0007585") # dim(GWAS.catalog.EFO.0007585) 42 38

#-------------------------------------------------------------------------------------------
# Run MR on CCPD and ICC cannabis initiation
#-------------------------------------------------------------------------------------------
library("TwoSampleMR",lib.loc="/mnt/backedup/home/lunC/R/x86_64-pc-linux-gnu-library/3.4")

CI.UKB.LDClumped.CCPD.MR.selected <- TwoSampleMR::mr(harmonised.CI.UKB.LDClumped.CCPD) %>% 
  filter(method %in% c("MR Egger","Weighted median","Inverse variance weighted"))    

CI.UKB.GCTACOJO.CCPD.MR.selected <- TwoSampleMR::mr(harmonised.CI.UKB.GCTACOJO.CCPD) %>% 
  filter(method %in% c("MR Egger","Weighted median","Inverse variance weighted"))    

#--------------------------------------------------------------------------------------------------
# Plot SNP effect on exposure against SNP effect on outcome by tweaking source code of TwoSampleMR::mr_scatter_plot()
#--------------------------------------------------------------------------------------------------
# Tweak the function TwoSampleMR::mr_scatter_plot() for a single plot

mr_scatter_plot_modified <-  function (mr_results
                                       ,dat
                                       ,colors
                                       ,exposure.label
                                       ,outcome.label){
  # mr_results: name of object created by TwoSampleMR::mr_()
  # dat: name of harmonised data created by TwoSampleMR::harmonise_data()
  # colors: User-supplied names of colors, one color per regression lines in the single plot
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
      ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure,y = beta.outcome)) + 
        ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome)
                           ,colour = "grey", width = 0) + 
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure)
                              ,colour = "grey", height = 0) + 
        ggplot2::geom_point(ggplot2::aes(text = paste("SNP:",SNP))) + 
        ggplot2::geom_abline(data = mrres
                             , ggplot2::aes(intercept = a, slope = b, colour = method)
                             , show.legend = TRUE) + 
        ggplot2::scale_colour_manual(values = colors) + 
        ggplot2::labs(colour = "Mendelian randomisation analyses with multiple SNPs"
                      , x = paste("SNP effect on", exposure.label) 
                      , y = paste("SNP effect on", outcome.label)) + 
        ggplot2::theme(legend.position = "top"
                       ,legend.direction = "vertical"
                       ,axis.text=element_text(size=14)) + 
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
      })
  mrres
}  

# Run the function
colors <- c("blue","black","orange")

plot.CI.UKB.LDClumped.CCPD <- mr_scatter_plot_modified(mr_results=CI.UKB.LDClumped.CCPD.MR.selected
                                      ,dat=harmonised.CI.UKB.LDClumped.CCPD
                                      ,colors=colors
                                      ,exposure.label="cups of coffee per day"
                                      ,outcome.label="cannabis initiation")

plot.CI.UKB.GCTACOJO.CCPD <- mr_scatter_plot_modified(mr_results=CI.UKB.GCTACOJO.CCPD.MR.selected
                                                       ,dat=harmonised.CI.UKB.GCTACOJO.CCPD
                                                       ,colors=colors
                                                       ,exposure.label="cups of coffee per day"
                                                       ,outcome.label="cannabis initiation")

## Export the scatter plots
ggplot2::ggsave(plot.CI.UKB.LDClumped.CCPD[[1]]
                , file=paste0(locPlots_MR,"zfig001-01_MR-SNP-effect_exposure-UKB-CCPD-LDClumped_outcome-ICC-cannabis-initiation.png")
                , width=7, height=7)

ggplot2::ggsave(plot.CI.UKB.GCTACOJO.CCPD[[1]]
                , file=paste0(locPlots_MR,"zfig001-02_MR-SNP-effect_exposure-UKB-CCPD-GCTACOJO_outcome-ICC-cannabis-initiation.png")
                , width=7, height=7)

#-------------------------------------------------------------------------------------------
# Sensitivity analyses- Heterogeneity statistics
## Hwo to interpret the result?
#-------------------------------------------------------------------------------------------
TwoSampleMR::mr_heterogeneity(harmonised.CI.UKB.LDClumped.CCPD)
TwoSampleMR::mr_heterogeneity(harmonised.CI.UKB.GCTACOJO.CCPD)

#-------------------------------------------------------------------------------------------
# Merge harmonised data and GWAS catalog download for trait coffee consumption ()
#-------------------------------------------------------------------------------------------
columns.keep <- c("SNPS"
                  ,"STRONGEST.SNP.RISK.ALLELE"
                  ,"MAPPED_TRAIT"
                  ,"REPORTED.GENE.S."
                  ,"MAPPED_GENE"
                  ,"UPSTREAM_GENE_ID")

GWAS.catalog.EFO.0004330.small <- GWAS.catalog.EFO.0004330 %>% select_(.dots=columns.keep)
GWAS.catalog.EFO.0007585.small <- GWAS.catalog.EFO.0007585 %>% select_(.dots=columns.keep)

harmoni.CI.UKB.LD.CCPD.catalog.EFO.0004330.short <- dplyr::left_join(harmonised.CI.UKB.LDClumped.CCPD
                                                               ,GWAS.catalog.EFO.0004330.small
                                                               ,by=c("SNP"="SNPS")) # dim(harmoni.CI.UKB.LD.CCPD.catalog.EFO.0004330.short) 19 35

harmoni.CI.UKB.LD.CCPD.catalog.EFO.0004330.full <- dplyr::left_join(harmonised.CI.UKB.LDClumped.CCPD
                                                                     ,GWAS.catalog.EFO.0004330
                                                                     ,by=c("SNP"="SNPS")) # dim(harmoni.CI.UKB.LD.CCPD.catalog.EFO.0004330.full) 19 67

harmoni.CI.UKB.COJO.CCPD.catalog.EFO.0004330.short <- dplyr::left_join(harmonised.CI.UKB.GCTACOJO.CCPD
                                                               ,GWAS.catalog.EFO.0004330.small
                                                               ,by=c("SNP"="SNPS")) # dim(harmoni.CI.UKB.COJO.CCPD.catalog.EFO.0004330.short) 19 67

harmoni.CI.UKB.COJO.CCPD.catalog.EFO.0004330.full <- dplyr::left_join(harmonised.CI.UKB.GCTACOJO.CCPD
                                                                       ,GWAS.catalog.EFO.0004330
                                                                       ,by=c("SNP"="SNPS")) # dim(harmoni.CI.UKB.COJO.CCPD.catalog.EFO.0004330.full) 19 67

#-------------------------------------------------------------------------------------#
#-----------------This is the end of this program-------------------------------------#
#-------------------------------------------------------------------------------------#
