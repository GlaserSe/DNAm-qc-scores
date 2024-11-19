#################################################################################
### Figure 5: Overlap differentially methylated CpGs (package: Upset)
#################################################################################

# set working directory
setwd("C:/Users/...")


# load packages
library(rio)
library(data.table)
library(dplyr)
library(ggplot2)
library(writexl)
library(UpSetR)


######################################################
# Import and prepare data
######################################################

# Import data frames for each differential DNA methylation analysis


# Create a list of all those dataframes
array_wo <- list("R1" = rotation_1, 
                 "R2" = rotation_2,
                 "R3" = rotation_3,
                 "R4" = rotation_4,
                 "R5" = rotation_5,
                 "R6" = rotation_6,
                 "R7" = rotation_7,
                 "R8" = rotation_8,
                 "R9" = rotation_9,
                 "R10" = rotation_10)



######################################################
# Overlap
######################################################

# Draw and save
pdf(file = "Figure 5_wo normalization control.pdf",
    width = 12,
    height = 7)
upset(fromList(array_wo), 
      nsets = 10,
      order.by = "freq",
      line.size = 1,
      mainbar.y.label = "CpG Intersections",
      sets.x.label = "Differentially methylated CpGs",
      text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 1.5),
      
      query.legend = "top",
      queries = list(list(query = intersects, params = list("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10"), 
                          color = "orange", 
                          active = T, 
                          query.name = "Included in all rotations"))
)
dev.off()