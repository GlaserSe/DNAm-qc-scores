###########################################################################################
# Quality assessment of DNA methylation data using Illumina Beadchip Arrays: BIN-Score
###########################################################################################
# last update, 30.01.2024


### Description of the function bin.score()

## Output
# Folder "BIN_Score" will be generated including a table comprising BIN-Scores of all samples

## Function arguments
# input: A Mset containing methylated and unmethylated signals (preferably generated with the minfi package). Reference samples must have "Control" in their name.
# array_type: Choose "EPIC" or "450k"


bin.score <- function(input, array_type){
  
  # Installing and loading required packages
  if (!requireNamespace("conumee")) install.packages("conumee")
  library(conumee)
  
  # Create folder for results
  dir.create("BIN_Score")
  
  results_CNV <- fread(paste0("Bin positions_", array_type, ".txt"), data.table=F)
  
  # Create annotation object
  data("exclude_regions")
  data("detail_regions")
  anno <- CNV.create_anno(array_type = array_type, 
                          exclude_regions = exclude_regions, 
                          detail_regions = detail_regions,
                          chrXY = TRUE)
  
  # Combine intensity values
  tumor.data    <- CNV.load(input)
  control.data  <- grep("Control", names(tumor.data))
  
  Mset <- mapToGenome(input)
  anno@probes<-subsetByOverlaps(anno@probes, granges(Mset))
  
  # Calculate the CNV-Score for each sample
  for(i in 1:ncol(tumor.data@intensity)) {
    
    # For each sample a CNV-Analysis will be performed
    sample  <- names(tumor.data@intensity)[i]
    x       <- CNV.fit(tumor.data[sample], tumor.data[control.data], anno)
    x       <- CNV.segment(CNV.detail(CNV.bin(x)))
    
    # Generation of a bin-list and alphabetical order of the bin-positions (chr)
    bin           <- CNV.write(x, what = "bins")
    binlist       <- fread(paste0("Bin positions_", array_type, ".txt"), data.table = F)
    binlist       <- cbind(binlist, bin)
    binlist       <- binlist[order(binlist$names),]
    bin           <- as.data.frame(bin[,5])
    
    # Generation of a segment-list  and alphabetical order of the seg-positions (chr)
    segment       <- CNV.write(x, what="segments")
    segmentlist   <- segment[,c(2,5,9)]                                       
    segmentlist   <- segmentlist[order(segmentlist$chrom),]
    
    # Generation of a vector including the corresponding number of segment median values for the bins within each segment
    anzahl_bins         <- as.vector(segmentlist$num.mark)    
    segmedian           <- as.vector(segmentlist$seg.median)  
    segmedian_all_bins  <- list()                             
    d                   <- 1                                  
    
    for(segment in 1:length(anzahl_bins)){
      
      for(anzahlbin in 1:anzahl_bins[segment]){
        wert                    <- segmedian[segment]
        segmedian_all_bins[[d]] <- wert
        d                       <- d +1
      }
    }
    
    # Combine bin-list with median values of segments
    segmedian_all_bins            <- as.data.frame(segmedian_all_bins)
    segmedian_all_bins            <- t(segmedian_all_bins)
    bin_median                    <- cbind(bin, segmedian_all_bins) 
    
    # Calculation difference of bins to segment line
    bin_median$diff     <-  bin_median$`bin[, 5]` - bin_median$segmedian_all_bins
    results_CNV             <- cbind(results_CNV, bin_median$diff)
    
  }
  
  # Adapt column names
  sample_names          <- as.vector(targets$Sample_Name)
  colnames(results_CNV) <- c("Bin positions", "Start", "End", "Width", sample_names)           
  
  # Calculation BIN-Score
  data        <- results_CNV[,5:ncol(results_CNV)]
  data_abs    <- abs(data)
  data_median <- apply(data_abs, 2, FUN=median, na.rm=TRUE)
  
  # Generation table including samples and corresponding BIN-Scores
  CNV                   <- as.data.frame(data_median)
  CNVscore              <- cbind(sample_names, CNV)
  colnames(CNVscore)    <- c("Sample", "BIN-Score")
  write.xlsx(CNVscore, "./BIN_Score/BIN Score.xlsx", colNames = T, rowNames=F)  
}


###  Example
bin.score(input = Mset, array_type = "450k")
