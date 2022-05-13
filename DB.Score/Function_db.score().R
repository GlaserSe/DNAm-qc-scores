###########################################################################################
# Quality assessment of DNA methylation data from Illumina Beadchip Arrays: DB-Score
###########################################################################################
# last update, 13.05.2022


### Description of the function db.score()

## Output
# Folder "DB_Score" will be generated including a table comprising DB-Scores of all samples

## Function arguments
# input: Name of the table with the TargetID and beta values (samples are listet per column)

db.score <- function(input){
  
  # Installing and loading required packages
  if (!requireNamespace("openxlsx")) install.packages("openxlsx")
  
  # Create folder for results
  dir.create("DB_Score")
  
  data_DBScore               <- input[,-1]
  DB                         <- list()
  hypo_list                  <- list()
  hyper_list                 <- list()
  mid_list                   <- list()
  
  for(Probe in 1:ncol(data_DBScore)) {
    
    # For each sample the number of CpGs <0.3, mid, >0.7 will be calculated
    hypo    <- nrow(subset(data_DBScore, data_DBScore[,Probe] < 0.3))
    hyper   <- nrow(subset(data_DBScore, data_DBScore[,Probe] > 0.7))
    mid     <- nrow(subset(data_DBScore, data_DBScore[,Probe] >= 0.3 & data_DBScore[,Probe] <= 0.7))
    
    # For each sample the QC-Score will be calculated and saved 
    score                <- round(mid/(hypo+hyper), 3)
    DB[[Probe]]          <- score
    hypo_list[[Probe]]   <- hypo
    hyper_list[[Probe]]  <- hyper
    mid_list[[Probe]]    <- mid
  }
  
  # Generating and saving a table including the DB-Scores of all samples
  DB                  <- t(as.data.frame(DB))
  hypo_list           <- t(as.data.frame(hypo_list))
  hyper_list          <- t(as.data.frame(hyper_list))
  mid_list            <- t(as.data.frame(mid_list))
  QC_table            <- as.data.frame(colnames(data_DBScore))                                                                 
  QC_table            <- cbind(QC_table, hypo_list, mid_list, hyper_list, DB)
  colnames(QC_table)  <- c("Sample", "< 0.3", "mid", "> 0.7", "DB-Score")
  write.xlsx(QC_table, "./DB_Score/DBScore.xlsx", colNames = T, rowNames=F) 
}


###  Example
db.score(input = data)
