###########################################################################################
# Quality assessment of DNA methylation data from Illumina Beadchip Arrays: CNV-Score
###########################################################################################
# last update, 13.05.2022


### Description of the function cnv.score()

## Output
# Folder "CNV_Score" will be generated including a table comprising CNV-Scores of all samples

## Function arguments
# input: A Mset containing methylated and unmethylated signals (preferably generated with the minfi package)
# array_type: Choose "EPIC" or "450k"


cnv.score <- function(input, array_type){
  
  # Installing and loading required packages
  if (!requireNamespace("conumee")) install.packages("conumee")
  library(conumee)
  
  # Create folder for results
  dir.create("CNV_Score")
  
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
    
    # Vektor generieren mit der entsprechende Anzahl der Seg.Median-Werte für die Bins in jedem Segment den Bins 
    # zugeordnet werden können 
    anzahl_bins         <- as.vector(segmentlist$num.mark)    # Vektor mit Anzahl der Bins in den einzelnen Segmenten
    segmedian           <- as.vector(segmentlist$seg.median)  # Vektor mit Segment-Mediane
    segmedian_all_bins  <- list()                             # Leere Liste, in der die Mediane hinzugefügt werden
    d                   <- 1                                  # Ermöglicht hinzufügen der Mediane (verhindert ersetzen)
    
    for(segment in 1:length(anzahl_bins)){
      
      for(anzahlbin in 1:anzahl_bins[segment]){
        wert                    <- segmedian[segment]
        segmedian_all_bins[[d]] <- wert
        d                       <- d +1
      }
    }
    
    # Verbinden der Probe-Binliste mit Mediane der Segmente
    segmedian_all_bins            <- as.data.frame(segmedian_all_bins)
    segmedian_all_bins            <- t(segmedian_all_bins)
    bin_median                    <- cbind(bin, segmedian_all_bins) 
    
    # Differenz der Bins zur blauen Segment-Linie berechnen
    bin_median$diff     <-  bin_median$`bin[, 5]` - bin_median$segmedian_all_bins
    results_CNV             <- cbind(results_CNV, bin_median$diff)
    
  }
  
  # Den Spalten die Proben-Namen zuweisen
  sample_names          <- as.vector(targets$Sample_Name)
  colnames(results_CNV) <- c("Bin positions", "Start", "End", "Width", sample_names)           
  
  # CNV-Score berechnen (Median)
  data        <- results_CNV[,5:ncol(results_CNV)]
  data_abs    <- abs(data)
  data_median <- apply(data_abs, 2, FUN=median, na.rm=TRUE)
  
  # Tabelle mit Proben + CNV-Score generieren
  CNV                   <- as.data.frame(data_median)
  CNVscore              <- cbind(sample_names, CNV)
  colnames(CNVscore)    <- c("Sample", "CNV-Score")
  write.xlsx(CNVscore, "./CNV_Score/CNV Score.xlsx", colNames = T, rowNames=F)  
}


###  Example
cnv.score(input = Mset, array_type = "450k")
