#################################################################################
### Figure 4: WGBS data
#################################################################################

# set working directory
setwd("C:/Users/...")


# load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(officer)
library(knitr)
library(tidyverse)



######################################################
# Import and prepare data
######################################################

# Import all data frames from WGBS and 450k data (containing beta-values)
# Sort by TargetID/rows
# Sort by Samples/columns
# Create a data frame for each sample
# Bring WGBS and 450k data to the sample TargetIDs



######################################################
# Figure 4A: Upper panel (beta values)
######################################################

# List of all data frames
df_list_diff <- list(
  "X4112512_cryo_mit"  = X4112512_cryo_mit, 
  "X4112512_cryo_ohne" = X4112512_cryo_ohne,
  "X4112512_FFPE_mit"  = X4112512_FFPE_mit, 
  "X4112512_FFPE_ohne" = X4112512_FFPE_ohne,
  
  "X4119027_cryo_mit"  = X4119027_cryo_mit, 
  "X4119027_cryo_ohne" = X4119027_cryo_ohne,
  "X4119027_FFPE_mit"  = X4119027_FFPE_mit, 
  "X4119027_FFPE_ohne" = X4119027_FFPE_ohne,
  
  "X4177434_cryo_mit"  = X4177434_cryo_mit, 
  "X4177434_cryo_ohne" = X4177434_cryo_ohne,
  "X4177434_FFPE_mit"  = X4177434_FFPE_mit, 
  "X4177434_FFPE_ohne" = X4177434_FFPE_ohne,
  
  "X4189998_cryo_mit"  = X4189998_cryo_mit, 
  "X4189998_cryo_ohne" = X4189998_cryo_ohne,
  "X4189998_FFPE_mit"  = X4189998_FFPE_mit, 
  "X4189998_FFPE_ohne" = X4189998_FFPE_ohne,
  
  "X4193278_cryo_mit"  = X4193278_cryo_mit, 
  "X4193278_cryo_ohne" = X4193278_cryo_ohne,
  "X4193278_FFPE_mit"  = X4193278_FFPE_mit, 
  "X4193278_FFPE_ohne" = X4193278_FFPE_ohne,
  
  "X4182393_cryo_mit"  = X4182393_cryo_mit, 
  "X4182393_cryo_ohne" = X4182393_cryo_ohne,
  "X4182393_FFPE_mit"  = X4182393_FFPE_mit, 
  "X4182393_FFPE_ohne" = X4182393_FFPE_ohne)

samples_wgbs <- c("4112512_WGBS", "4112512_WGBS", "4112512_WGBS", "4112512_WGBS",
                  "4119027_WGBS", "4119027_WGBS", "4119027_WGBS", "4119027_WGBS",
                  "4177434_WGBS", "4177434_WGBS", "4177434_WGBS", "4177434_WGBS",
                  "4193278_WGBS", "4193278_WGBS", "4193278_WGBS", "4193278_WGBS",
                  "4182393_WGBS", "4182393_WGBS", "4182393_WGBS", "4182393_WGBS",
                  "4189998_WGBS", "4189998_WGBS", "4189998_WGBS", "4189998_WGBS")


# Loop to create density plots
for (i in 1:24) {
  
  # Select wgbs and array data, combine in one table
  wgbs_array  <- inner_join(wgbs[,c("TargetID", samples_wgbs[i])], df_list[[i]], by = "TargetID")
  wgbs_array  <- wgbs_array[,-c(1)]
  
  # Adapt table for plot generation
  colnames(wgbs_array)  <- c("WGBS", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
  wgbs_array            <- wgbs_array %>% mutate_if(is.character, as.numeric)
  wgbs_array_long       <- gather(wgbs_array, Group, Value, "WGBS":"R10")
  wgbs_array_long$Group <- factor(wgbs_array_long$Group,
                                  levels = c("WGBS", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10"))
  
  # Plot
  plot_betaval <- ggplot(wgbs_array_long, aes(Value, fill =Group, color = Group)) + 
    geom_density(alpha = 0.3) +
    scale_color_manual(values = c("#67000d", "#737373", "#737373", "#737373", "#737373", "#737373", "#737373", "#737373", "#737373", "#737373", "#737373"))+
    scale_fill_manual(values = c("#b30000", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd"))+
    ggtitle(names(df_list)[[i]]) +
    scale_x_continuous(name = "Beta value") +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_blank(), 
          axis.line = element_line(color="black"), 
          axis.line.x = element_line(color="black"),
          legend.title = element_text(face = "bold"))
  
  # Save
  ggsave(filename = paste0("Figure 4A_Beta value_WGBS_vs_", names(df_list_diff)[[i]], ".pdf"),
         plot = plot_betaval,
         width = 7,
         height = 4)
}



######################################################
# Figure 4A: Lower panel (difference)
######################################################

# List of all dataframes
df_list_diff <- list(
  "X4112512_cryo_mit"  = X4112512_cryo_mit, 
  "X4112512_cryo_ohne" = X4112512_cryo_ohne,
  "X4119027_cryo_mit"  = X4119027_cryo_mit, 
  "X4119027_cryo_ohne" = X4119027_cryo_ohne,
  "X4177434_cryo_mit"  = X4177434_cryo_mit, 
  "X4177434_cryo_ohne" = X4177434_cryo_ohne,
  "X4189998_cryo_mit"  = X4189998_cryo_mit, 
  "X4189998_cryo_ohne" = X4189998_cryo_ohne,
  "X4193278_cryo_mit"  = X4193278_cryo_mit, 
  "X4193278_cryo_ohne" = X4193278_cryo_ohne,
  "X4182393_cryo_mit"  = X4182393_cryo_mit, 
  "X4182393_cryo_ohne" = X4182393_cryo_ohne)

samples_wgbs_diff <- c("4112512_WGBS", "4112512_WGBS", 
                       "4119027_WGBS", "4119027_WGBS", 
                       "4177434_WGBS", "4177434_WGBS", 
                       "4193278_WGBS", "4193278_WGBS", 
                       "4182393_WGBS", "4182393_WGBS", 
                       "4189998_WGBS", "4189998_WGBS")


# Create final table with absolute mean difference
df_absmeandiff <- data.frame()
for (i in 1:12) {
  
  # select wgbs and array data, combine in one table
  wgbs_array            <- inner_join(wgbs[,c("TargetID", samples_wgbs_diff[i])], df_list_diff[[i]], by = "TargetID")
  wgbs_array            <- wgbs_array[,-c(1)]
  colnames(wgbs_array)  <- c("WGBS", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10")
  wgbs_array            <- wgbs_array %>% mutate_if(is.character, as.numeric)
  
  # calculate differences
  wgbs_array <- mutate(wgbs_array, Minus_R1 = abs(WGBS - R1))
  wgbs_array <- mutate(wgbs_array, Minus_R2 = abs(WGBS - R2))
  wgbs_array <- mutate(wgbs_array, Minus_R3 = abs(WGBS - R3))
  wgbs_array <- mutate(wgbs_array, Minus_R4 = abs(WGBS - R4))
  wgbs_array <- mutate(wgbs_array, Minus_R5 = abs(WGBS - R5))
  wgbs_array <- mutate(wgbs_array, Minus_R6 = abs(WGBS - R6))
  wgbs_array <- mutate(wgbs_array, Minus_R7 = abs(WGBS - R7))
  wgbs_array <- mutate(wgbs_array, Minus_R8 = abs(WGBS - R8))
  wgbs_array <- mutate(wgbs_array, Minus_R9 = abs(WGBS - R9))
  wgbs_array <- mutate(wgbs_array, Minus_R10 = abs(WGBS - R10))
  wgbs_array <- wgbs_array[,12:21]
  
  # adapt table for plot
  colnames(wgbs_array)
  wgbs_array_long       <- gather(wgbs_array, Group, Value, "Minus_R1":"Minus_R10")
  
  
  # plot density
  plot_diff_dens <- ggplot(wgbs_array_long, aes(Value, fill =Group, color = Group)) + 
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a", "#1a1a1a"))+
    scale_fill_manual(values = c("#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#bdbdbd"))+
    ggtitle(names(df_list_diff)[[i]]) +
    scale_x_continuous(name = "Beta value") +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_blank(), 
          axis.line = element_line(color="black"), 
          axis.line.x = element_line(color="black"),
          legend.title = element_text(face = "bold"))
  
  # save
  ggsave(filename = paste0("Figure 4A_Difference_WGBS_vs_", names(df_list_diff)[[i]], ".pdf"),
         plot = plot_diff_dens,
         width = 7,
         height = 4)
  
  
  # calculate absolute mean difference per sample
  wgbs_array_mean <- wgbs_array %>% 
    colMeans(na.rm = T) %>% 
    as.data.frame() %>% 
    setDT(keep.rownames = TRUE)
  colnames(wgbs_array_mean) <- c("Group", "Abs.mean")
  wgbs_array_mean$Sample <- names(df_list_diff)[[i]]
  
  # merge with final df
  df_absmeandiff <- rbind(df_absmeandiff, wgbs_array_mean)
}


######################################################
# Figure 4B: Scatter plot absolute mean difference
######################################################

# Prepare table for plot
df_absmeandiff <- as.data.frame(df_absmeandiff)
df_absmeandiff[c('Sample', 'Material', 'Method')] <- str_split_fixed(df_absmeandiff$Sample, '_', 3)
df_absmeandiff$Sample <- sub('.', '', df_absmeandiff$Sample)
df_absmeandiff <- df_absmeandiff  %>%
  mutate(across('Method', str_replace, 'mit', 'with normalization samples'))
df_absmeandiff <- df_absmeandiff  %>%
  mutate(across('Material', str_replace, 'cryo', 'cryo-preserved'))
df_absmeandiff <- df_absmeandiff  %>%
  mutate(across('Method', str_replace, 'ohne', 'without normalization samples'))
df_absmeandiff$Method = factor(df_absmeandiff$Method, levels=c('without normalization samples','with normalization samples'))


# plot
plot_absdiff <- ggplot(df_absmeandiff, aes(x=Sample, y=Abs.mean)) +
  geom_jitter(width = 0.2) +
  xlab("") +
  ylab("absolute mean difference") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12)) +
  facet_wrap(~ Method)

# save
ggsave(filename = "Figure 4B_Scatter plot difference.pdf",
       plot = plot_absdiff,
       width = 6,
       height = 4)



######################################################
# Plots for Figure 4A were further adapted with the program Inkscape
######################################################
