#################################################################################
### Figure 6: Application to different sample types
#################################################################################

# set working directory
setwd("C:/Users/...")


# load packages
library(rio)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggbreak)
library(ggforce)
library(plotrix)
library(data.table)
library(dplyr)


######################################################
# Example Figure 6A
######################################################

# Import (dataframe containing values for BIN- and DB-Score)
data <- import("Plot_BIN vs DB.csv")


# Plot
ggplot(data = data, aes(x= `QC.Score`, y = `CNV.Score`, color = `Quality`, label = `Sample`)) +
  geom_point(size=4) +
  scale_color_manual(values=c("#a50f15","#525252","#7fbc41", "#fe9929", "red")) + 
  xlab("DB-Score") + 
  ylab("BIN-Score") +
  geom_vline(xintercept = 1, colour = "red", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.25, colour = "red", linewidth = 1, linetype = "dashed") +
  theme_classic(20) +
  theme(axis.line.x   = element_line(colour="black"), 
        axis.line.y   = element_line(colour="black"), 
        panel.spacing = unit(1, "lines"),
        axis.title.x  = element_text(size =16),
        axis.title.y  = element_text(size =16),
        legend.position = "right",
        legend.title = element_text(face="bold")) +
  facet_zoom(xlim = c(0, 2),
             ylim = c(0, 0.5),
             horizontal = F,
             zoom.size = 4,
             shrink = F)

# Save
ggsave(filename = "Figure 6A_BIN vs DB.pdf",
       width = 20,
       height = 15)



######################################################
# Example Figure 6B
######################################################


# Plot
ggplot(data, aes(x=`Array_4115001`, y=`WGBS_4115001`)) + 
  geom_point(alpha= 0.5, color = "darkgrey") +
  stat_density_2d(geom = "polygon", contour = T,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 7)+
  scale_fill_distiller(palette = "Reds") +
  scale_x_continuous(name = "450k", 
                     limits = c(0, 1), 
                     breaks = seq(0, 1, 0.5)) + 
  scale_y_continuous(name = "WGBS", 
                     limits = c(0,1), 
                     breaks = seq(0, 1, 0.5)) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_blank(), 
        axis.line = element_line(color="black"), 
        axis.text = element_text(size = 20, color = "black")) +
  stat_cor(method="pearson", size = 7, label.sep = "\n", label.x = 0, label.y = 0.95) 


# Save
ggsave(filename = "Figure 6B_Scatter density.png",
       dpi = 300,
       height = 5, 
       width = 10) 
