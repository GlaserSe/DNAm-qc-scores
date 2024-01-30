###########################################################################################
# Quality assessment of DNA methylation data using Illumina Beadchip Arrays: CM-Scores
###########################################################################################
# last update, 30.01.2024


### Description of the function cm.score()

## Output
# Folder "CM_Score" will be generated including a table comprising CM-Scores of all samples

## Function arguments
# input: Name of the table with the TargetID and beta values (samples are listed per column)

cm.score <- function(input){
  
  # Installing and loading required packages
  if (!requireNamespace("data.table")) install.packages("data.table")
  library(data.table)
  if (!requireNamespace("dplyr")) install.packages("dplyr")
  library(dplyr)
  
  # Create results folder and load stable loci
  dir.create("CM_Score")
  stableloci   <- fread("Stable_loci.txt")
  sample_names <- as.vector(colnames(input[,2:ncol(input)]))
  
  # Filter within your input table for the stableloci
  data_CM             <- inner_join(stableloci,input, by = "TargetID") 
  data_CM             <- as.data.frame(data_CM)
  
  # CM Score calculation
  CMvalue             <- list()
  CMhigh              <- list()
  CMlow               <- list()
  
  for (sample in 1:ncol(data_CM)) {
    high_good <- nrow(subset(data_CM, data_CM$Type == "high" & data_CM[,sample] > 0.9))
    low_good  <- nrow(subset(data_CM, data_CM$Type == "low" & data_CM[,sample] < 0.1))
    
    CMscore_high   <- round(100-((high_good/279)*100),2) 
    CMscore_low    <- round(100-((low_good/313)*100),2)
    values         <- c(CMscore_high, CMscore_low)
    CMscore        <- mean(values)
    
    CMhigh[[sample]]  <- CMscore_high
    CMlow[[sample]]   <- CMscore_low
    CMvalue[[sample]] <- CMscore
  }
  
  # Scores als Dataframes speichern (Cave: Erste zwei spalten waren TargetID und Type --> l?schen)
  CMhigh          <- t(as.data.frame(CMhigh))
  CMhigh          <- as.data.frame(CMhigh[3:ncol(data_CM)])
  CMlow           <- t(as.data.frame(CMlow))
  CMlow           <- as.data.frame(CMlow[3:ncol(data_CM)])
  
  # Save Scores in a table
  data_hilo             <- cbind(sample_names, CMhigh, CMlow)
  colnames(data_hilo)   <- c("Sample", "CM.Score.high", "CM.Score.low")
  
  
  # Calculate HL.Score.diff
  data_hilo$CM.Score.diff        <- abs(data_hilo$CM.Score.high - data_hilo$CM.Score.low)
  
  # Save results tabel in the above created folder
  write.table(data_hilo, "./CM_Score/CMScores.csv", row.names = F, sep = ",")
}


###  Example
cm.score(input=data)