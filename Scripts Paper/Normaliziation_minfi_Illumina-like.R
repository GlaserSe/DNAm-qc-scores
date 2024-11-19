#############################################################################################
#############################################################################################
# Methylierungsanalyse mittels Minfi
#############################################################################################
#############################################################################################

#+++++++++++++++++++++++++++++++++++++++++++++++++++
path_idat                   <- "X:/folder_with_all_idats"
path_results                <- "C:/Users/results"
array_type                  <- "450k"   
#+++++++++++++++++++++++++++++++++++++++++++++++++++

# set working directory
setwd(path_results)


# load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(minfi)
library(conumee2.0)
library(rlist)
library(tidyr)
library(magrittr)
library(writexl)
library(officer)
library(ComplexHeatmap)
library(circlize)
library(rio)


# import files
rotations   <- import("Describing_alll_rankings.xlsx")
samplesheet <- import("SampleSheet.xlsx")


# Loop to perform normalization through all sample sheet rankings
for (i in 1:ncol(rotations)) {
  
  #############################################################################################
  # import data and normalize
  #############################################################################################
  
  # generate samplesheet
  sample_order        <- rotations[,i]
  targets             <- samplesheet[match(sample_order, samplesheet$Sample_Name),]
  targets$Basename    <- paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)                                                    
  RGSet               <- read.metharray.exp(base = path_idat, targets=targets, force = TRUE, extended = TRUE, recursive = TRUE)
  targets$ID          <- paste(targets$Sample_Name)                                                                   
  sampleNames(RGSet)  <- targets$ID
  
  # Illumina-like normalization
  Mset    <- preprocessIllumina(RGSet, normalize = "controls", bg.correct = FALSE)
  GMset   <- mapToGenome(Mset)
  methy   <- minfi::getMeth(GMset)
  unmethy <- minfi::getUnmeth(GMset)
  betas   <- methy / (methy +unmethy +100)
  betas   <- as.data.frame(betas)
  betas   <- betas[,-1]
  
  
  #############################################################################################
  # Detection p-value filtering  (< 0.01)
  #############################################################################################
  
  # data frame with detection p-values
  DetSet <- as.data.frame(detectionP(RGSet, type="m+u"))
  DetSet <- round(DetSet, 6)
  bval_table    <- betas
  pval_table    <- DetSet[,2:ncol(DetSet)]                                                                           
  setDT(pval_table, keep.rownames = TRUE)
  colnames(pval_table)[1] <- "TargetID"
  setDT(bval_table, keep.rownames = TRUE)
  colnames(bval_table)[1] <- "TargetID"
  
  # change to long format
  bval_table_long <- bval_table %>% 
    gather(Sample, B_value, -TargetID)
  pval_table_long <- pval_table %>% 
    gather(Sample, p_value, -TargetID)
  
  # merge beta-values and p-values based on sample name and TargetID
  bpval_merged <- merge(bval_table_long, pval_table_long, by = c("Sample", "TargetID"))
  
  # P-value filtering (< 0.01)
  pval_filtered <- bpval_merged %>% 
    mutate(B_val_new = ifelse(p_value > 0.01, NA, B_value))
  
  pval_filtered <- pval_filtered[,c(1,2,5)]

  # change to wide format
  pval_filtered_wide <- pval_filtered %>% 
    spread(Sample, B_val_new)
  betas_pdet <- pval_filtered_wide
  order_samples    <- targets$Sample_Name[-(1)]                                                         
  order_samples    <- c("TargetID", order_samples)
  betas_pdet_order <- betas_pdet[,order_samples]
  
  #############################################################################################
  # Save files
  #############################################################################################
  
  # Import array annotation
  info <- import("")
  
  # Merge beta-values and array annotation
  betas_pdet_infos   <- inner_join(info, betas_pdet_order, by = "TargetID") 
  
  # exclude X/Y-Loci 
  betas_pdet_infos_wo_X <- betas_pdet_infos[-grep("^X",betas_pdet_infos$CHR),]
  betas_pdet_infos_wo_XY <- betas_pdet_infos_wo_X[-grep("^Y",betas_pdet_infos_wo_X$CHR),]
  
  # save
  write.table(betas_pdet_infos_wo_XY, ".txt", row.names = F)
}


