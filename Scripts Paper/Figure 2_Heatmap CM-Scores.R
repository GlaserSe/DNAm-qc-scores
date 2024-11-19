##################################################################
# Figure 2: Heatmap CM-Scores
##################################################################

# Set working directory
setwd("C:/Users/...")


# Load packages
library(data.table)
library(circlize)
library(colorspace)
library(GetoptLong)
library(ComplexHeatmap)
library(rio)



######################################################
# Import and prepare data
######################################################

# Import data (data frame containing beta-values. columns: samples, rows: TargetIDs)
data   <- import("File.csv")
colnames(data)


# Filter data and create Matrix
hm_data            <- data
row.names(hm_data) <- hm_data[,1]
hm_data            <- hm_data[,-1]
hm_data_mat        <- as.matrix(hm_data)



######################################################
# Create Heatmap (package: ComplexHeatmap)
######################################################

# Import Annotation (Cave: Columns of matrix and annotation have to match!!)
anno   <- import("Annotation.csv")


# Column annotation
col_ha = HeatmapAnnotation("Legend" = anno$Type,
                           col = list("Legend" = c("normalization sample"     = "#662506",
                                                   "bad CM-Score difference"  = "#ec7014",
                                                   "bad CM-Score high and low"= "#fec44f")))

# Row annotation
row_anno # define row annotations
row_ha = rowAnnotation("Stable loci" = row_anno,
                       col = list("Stable loci" = c("high" = "black",
                                                    "low"  = "#bdbdbd")))


# Heatmap: normalization samples
ht1 <- Heatmap(hm_data_mat[,1:5],
               col = colorRamp2(c(0,0.5,1),
                                c("blue", "black", "yellow")), 
               name = "DNA methylation (beta-value)",
               show_row_names = FALSE,
               show_column_names = FALSE,
               left_annotation = row_ha,
               cluster_columns = FALSE,
               cluster_rows = F,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           legend_width = unit(5, "cm"))) 
draw(ht1,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE,
     padding = unit(c(10,10,10,10), "mm"))


# Heatmap: bad CM-Score difference
ht2 <- Heatmap(hm_data_mat[,6:8],
               col = colorRamp2(c(0,0.5,1),
                                c("blue", "black", "yellow")), 
               name = "DNA methylation (beta-value)",
               show_row_names = FALSE,
               show_column_names = FALSE,
               left_annotation = row_ha,
               cluster_columns = FALSE,
               cluster_rows = F,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           legend_width = unit(5, "cm"))) 
draw(ht2,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE,
     padding = unit(c(10,10,10,10), "mm"))


# Heatmap: bad CM-Score high and low
ht3 <- Heatmap(hm_data_mat[,9:13],
               col = colorRamp2(c(0,0.5,1),
                                c("blue", "black", "yellow")), 
               name = "DNA methylation (beta-value)",
               show_row_names = FALSE,
               show_column_names = FALSE,
               left_annotation = row_ha,
               cluster_columns = FALSE,
               cluster_rows = F,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           legend_width = unit(5, "cm"))) 
draw(ht3,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE,
     padding = unit(c(10,10,10,10), "mm"))


# Draw all three heatmaps
ht_list = ht1 + ht2 + ht3 


# Save
pdf(file = "Figure 2_Heatmap CM Scores.pdf", 
    width = 10, 
    height = 10) 
draw(ht_list,
     annotation_legend_side = "bottom",
     heatmap_legend_side = "bottom",
     merge_legend = TRUE,
     padding = unit(c(10,10,10,10), "mm"))
dev.off()


######################################################
# Afterwards the figure was adapted with the program Inkscape
######################################################
