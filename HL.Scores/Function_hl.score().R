###########################################################################################
# Quality assessment of DNA methylation data from Illumina Beadchip Arrays: HL-Scores
###########################################################################################
# last update, 13.05.2022


### Description of the function hl.score()

## Output
# Folder "HL_Score" will be generated including a table comprising HL-Scores of all samples

## Function arguments
# input: Name of the table with the TargetID and beta values (samples are listed per column)

hl.score <- function(input){
  
  # Installing and loading required packages
  if (!requireNamespace("data.table")) install.packages("data.table")
  library(data.table)
  if (!requireNamespace("dplyr")) install.packages("dplyr")
  library(dplyr)
  
  # Create results folder and load stable loci
  dir.create("HL_Score")
  stableloci <- fread("Stable_loci.txt")
  sample_names <- as.vector(colnames(input[,2:ncol(input)]))
  
  # Filter within your input table for the stableloci
  data_HL             <- inner_join(stableloci,input, by = "TargetID") 
  data_HL             <- as.data.frame(data_HL)
  
  # HL Score calculation
  HLvalue             <- list()
  HLhigh              <- list()
  HLlow               <- list()
  
  for (sample in 1:ncol(data_HL)) {
    high_good <- nrow(subset(data_HL, data_HL$Type == "high" & data_HL[,sample] > 0.9))
    low_good  <- nrow(subset(data_HL, data_HL$Type == "low" & data_HL[,sample] < 0.1))
    
    HLscore_high   <- round(100-((high_good/279)*100),2) 
    HLscore_low    <- round(100-((low_good/313)*100),2)
    values         <- c(HLscore_high, HLscore_low)
    HLscore        <- mean(values)
    
    HLhigh[[sample]]  <- HLscore_high
    HLlow[[sample]]   <- HLscore_low
    HLvalue[[sample]] <- HLscore
  }
  
  # Scores als Dataframes speichern (Cave: Erste zwei spalten waren TargetID und Type --> löschen)
  HLhigh          <- t(as.data.frame(HLhigh))
  HLhigh          <- as.data.frame(HLhigh[3:ncol(data_HL)])
  HLlow           <- t(as.data.frame(HLlow))
  HLlow           <- as.data.frame(HLlow[3:ncol(data_HL)])
  
  # Save Scores in a table
  data_hilo             <- cbind(sample_names, HLhigh, HLlow)
  colnames(data_hilo)   <- c("Sample", "HL.Score.high", "HL.Score.low")
  
  
  # Calculate HL.Score.diff
  data_hilo$HL.Score.diff        <- abs(data_hilo$HL.Score.high - data_hilo$HL.Score.low)
  
  # Save results tabel in the above created folder
  write.table(data_hilo, "./HL_Score/HLScores.csv", row.names = F, sep = ",")
}


###  Example
hl.score(input=data)
