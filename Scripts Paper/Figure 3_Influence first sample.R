#################################################################################
### Figure 3: Influence of the first listed sample using Illumina-like normalization algorithm
#################################################################################

# Set working directory
setwd("C:/Users/...")


# load packages
library(rio)
library(tidyverse)
library(stringr)
library(tidyr)
library(matrixStats)
library(data.table)
library(GGally)
library(matrixStats)
library(stringr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggbreak)
library(ggforce)



######################################################
# Import and prepare data
######################################################

# Import all data frames (containing beta-values, normalized with and without normalization sample)
# Sort by TargetID/rows
# Sort by Samples/columns
# Create a data frame for each sample (e.g., X4112512_cryo_ohne, X4112512_cryo_mit, ..)


# Create a list with all dataframes
df_list      <- list(X4112512_cryo_ohne, X4112512_cryo_mit, 
                     X4112512_FFPE_ohne, X4112512_FFPE_mit,
                     X4119027_cryo_ohne, X4119027_cryo_mit,
                     X4119027_FFPE_ohne, X4119027_FFPE_mit,
                     X4177434_cryo_ohne, X4177434_cryo_mit, 
                     X4177434_FFPE_ohne, X4177434_FFPE_mit,
                     X4182393_cryo_ohne, X4182393_cryo_mit, 
                     X4182393_FFPE_ohne, X4182393_FFPE_mit,
                     X4189998_cryo_ohne, X4189998_cryo_mit, 
                     X4189998_FFPE_ohne, X4189998_FFPE_mit,
                     X4193278_cryo_ohne, X4193278_cryo_mit, 
                     X4193278_FFPE_ohne, X4193278_FFPE_mit)
names(df_list) <-  c("X4112512_cryo_ohne", "X4112512_cryo_mit", 
                     "X4112512_FFPE_ohne", "X4112512_FFPE_mit",
                     "X4119027_cryo_ohne", "X4119027_cryo_mit",
                     "X4119027_FFPE_ohne", "X4119027_FFPE_mit",
                     "X4177434_cryo_ohne", "X4177434_cryo_mit", 
                     "X4177434_FFPE_ohne", "X4177434_FFPE_mit",
                     "X4182393_cryo_ohne", "X4182393_cryo_mit", 
                     "X4182393_FFPE_ohne", "X4182393_FFPE_mit",
                     "X4189998_cryo_ohne", "X4189998_cryo_mit", 
                     "X4189998_FFPE_ohne", "X4189998_FFPE_mit",
                     "X4193278_cryo_ohne", "X4193278_cryo_mit", 
                     "X4193278_FFPE_ohne", "X4193278_FFPE_mit")



######################################################
# Figure 3A-B: Correlation matrix
######################################################


for (i in 1:length(df_list)) {
  
  # select sample and convert to matrix
  data           <- as.data.frame(df_list[[i]])
  colnames(data) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
  
  
  # save plots
  tiff(filename = paste0("Correlation_Pearson_", names(df_list)[i], ".tiff"),
       width = 1200,
       height = 1000)
  plot <- ggpairs(data,
                  upper = list(continuous = wrap(ggally_cor, method = "pearson", corSize=15))) +
    theme_bw() +
    theme(strip.text = element_text(size = 12))
  print(plot)
  dev.off()
}



######################################################
# Figure 3C: Plot Pearson correlation
######################################################


# Create table with all correlation coefficients
matrix           <- as.matrix(df_list[[1]])
colnames(matrix) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
cormat           <- round(cor(matrix, method = "spearman", use = "complete.obs"), 2)
melted_cormat    <- melt(cormat)
df_corr_val      <- melted_cormat[,c(1:2)]


for (i in 1:length(df_list)) {
  
  # select sample and convert to matrix
  matrix           <- as.matrix(df_list[[i]])
  colnames(matrix) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
  
  # correlation
  cormat        <- round(cor(matrix, method = "pearson", use = "complete.obs"), 2)
  melted_cormat <- melt(cormat)
  
  # Prepare table
  colnames(melted_cormat) <- c("Var1", "Var2", names(df_list[i]))
  
  # add to df_corr_val
  df_corr_val <- merge(df_corr_val, melted_cormat, by = c("Var1", "Var2"))
}


# Prepare dataframe for plot
df_plot_corr <- gather(df_corr_val, Sample, Value, "X4112512_cryo_ohne":"X4193278_FFPE_mit")
df_plot_corr[c('Sample', 'Fix.Method', 'Norm')] <- str_split_fixed(df_plot_corr$Sample, '_', 3)
df_plot_corr$Sample <- sub('.', '', df_plot_corr$Sample)
df_plot_corr <- df_plot_corr  %>%
  mutate(across('Norm', str_replace, 'mit', 'with normalization samples'))
df_plot_corr <- df_plot_corr  %>%
  mutate(across('Fix.Method', str_replace, 'cryo', 'cryo-preserved'))
df_plot_corr <- df_plot_corr  %>%
  mutate(across('Norm', str_replace, 'ohne', 'without normalization samples'))
df_plot_corr$Norm = factor(df_plot_corr$Norm, levels=c('without normalization samples','with normalization samples'))


# Box plot
ggplot(df_plot_corr, aes(x=Sample,y=Value)) + 
  geom_jitter(width = 0.2) +
  xlab("") +
  ylab("Pearson correlation coefficient R") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12)) +
  facet_grid(Fix.Method ~ Norm) 


# Save
ggsave(filename = "Figure 3C_Boxplot_Pearson.pdf",
       width = 8,
       height = 6)




######################################################
# Figure 3D: Plot Standard deviation
######################################################

# Create dataframe to save results (final df)
df_var_final    <- data.frame(R1_mit$TargetID)

for (i in 1:length(df_list)) {
  # calculate SD per CpG
  df_var           <- as.data.frame(df_list[i])
  colnames(df_var) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
  df_var           <- as.data.frame(apply(df_var, 1, sd, na.rm=T))
  colnames(df_var) <- names(df_list)[[i]]
  
  # merge with final df
  df_var_final <- cbind(df_var_final, df_var)
}


# Prepare dataframe for plot
colnames(df_var_final)[1] <- "TargetID"
df_var_final_long <- gather(df_var_final, Sample, SD, "X4112512_cryo_ohne":"X4193278_FFPE_mit")
df_var_final_long[,c("Sample", "Material", "Method")] <- str_split_fixed(df_var_final_long$Sample, "_", 3)
df_var_final_long$Sample <- sub('.', '', df_var_final_long$Sample)
df_var_final_long <- df_var_final_long  %>%
  mutate(across('Method', str_replace, 'mit', 'with normalization samples'))
df_var_final_long <- df_var_final_long  %>%
  mutate(across('Material', str_replace, 'cryo', 'cryo-preserved'))
df_var_final_long <- df_var_final_long  %>%
  mutate(across('Method', str_replace, 'ohne', 'without normalization samples'))
df_var_final_long$Method = factor(df_var_final_long$Method, levels=c('without normalization samples','with normalization samples'))


# Plot
ggplot(df_var_final_long, aes(x=Sample, y=SD)) +
  geom_boxplot() +
  xlab("") +
  ylab("standard deviation across rotations") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12)) +
  facet_grid(Material ~ Method) 


# Save
ggsave(filename = "Figure 3D_Boxplot_standard deviation.pdf",
       width = 8,
       height = 6)




######################################################
# Figure 3E-F: Bar plots differing CpGs
######################################################

### Filter number of differing CpGs
data_diff_01_final <- data.frame()

for (i in 1:length(df_list)) {
  
  data_diff             <- df_list[[i]]
  data_diff_01           <- rep(c(names(df_list[i])),10)
  data_diff_01           <- as.data.frame(data_diff_01)
  data_diff_01$Rotation  <- c("R1 minus R1", "R1 minus R2", "R1 minus R3", "R1 minus R4", "R1 minus R5", "R1 minus R6", "R1 minus R7", "R1 minus R8", "R1 minus R9", "R1 minus R10")
  data_diff_01$Numb.Diff <- "NA"
  
  for(r in 1:ncol(data_diff)){
    R1                <- as.data.frame(data_diff[,1])
    Rother            <- as.data.frame(data_diff[,r]) 
    diff           <- cbind(R1, Rother)
    colnames(diff) <- c("R1", "Diff")
    diff           <- transmute(diff, name = abs(R1-Diff))
    data_diff_01[r,3]      <- nrow(subset(diff, name > 0.1))
    table_save             <- subset(diff, name > 0.1)
    write.table(table_save, paste0("Tables_diff_0.1/Diff0.1_", names(df_list[i]), "_R1_vs_R", r, ".csv"), row.names = F, sep = ",")
  }
  data_diff_01_final <- rbind(data_diff_01_final, data_diff_01)
}


# split first column
data_diff_01_final[c('Sample', 'Material', 'Process')] <- str_split_fixed(data_diff_01_final$data_diff_01, '_', 3)
data_diff_01_final <- data_diff_01_final  %>%
  mutate(across('Process', str_replace, 'mit', 'with normalization sample'))
data_diff_01_final <- data_diff_01_final  %>%
  mutate(across('Process', str_replace, 'ohne', 'without normalization sample'))


### Plot
# New plot: Stacked bar plot
data_plot           <- data_diff_01_final
data_plot$Numb.Diff <- as.numeric(data_plot$Numb.Diff)
data_plot_wo        <- subset(data_plot, data_plot$Process == "without normalization sample")
data_wo_cryo        <- subset(data_plot_wo, data_plot_wo$Material == "cryo")
data_wo_ffpe        <- subset(data_plot_wo, data_plot_wo$Material == "FFPE")
data_wo_cryo_wo     <- subset(data_wo_cryo, data_wo_cryo$Number.CpGs > 0)
data_wo_ffpe_wo     <- subset(data_wo_ffpe, data_wo_ffpe$Number.CpGs > 0)


# Cryo
plot_cryo <- 
  ggplot(data_wo_cryo, aes(x=Rotation, y=Numb.Diff, fill=Sample)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("#c7e9c0", "#41ab5d", "#006d2c", "#9ecae1", "#4292c6", "#08519c")) +
  scale_y_continuous(limits = c(-1,280000)) +
  scale_x_discrete(limits = c("R1 minus R2", "R1 minus R3", "R1 minus R4", "R1 minus R5","R1 minus R6", "R1 minus R7", "R1 minus R8", "R1 minus R9", "R1 minus R10")) +
  labs(x="",
       y="Number of differing CpGs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold")) +
  facet_zoom(ylim = c(0, 10000))
print(plot_cryo)


# Save
ggsave(filename = "Figure 3D_Cryo.pdf",
       plot = plot_cryo,
       width = 9,
       height = 4)


# FFPE
plot_ffpe <- 
  ggplot(data_wo_ffpe, aes(x=Rotation, y=Numb.Diff, fill=Sample)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("#c7e9c0", "#41ab5d", "#006d2c", "#9ecae1", "#4292c6", "#08519c")) +
  scale_y_continuous(limits = c(-1,280000)) +
  scale_x_discrete(limits = c("R1 minus R2", "R1 minus R3", "R1 minus R4", "R1 minus R5","R1 minus R6", "R1 minus R7", "R1 minus R8", "R1 minus R9", "R1 minus R10")) +
  labs(x="",
       y="Number of differing CpGs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold")) +
  facet_zoom(ylim = c(0, 10000))
print(plot_ffpe)


# Save
ggsave(filename = "Figure 3D_FFPE.pdf",
       plot = plot_ffpe,
       width = 9,
       height = 4)
