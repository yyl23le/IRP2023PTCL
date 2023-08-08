# Authors: Friday
# Modified Date: 30-05-2023
# Version: 1.0

# The following script analyse the fragmentation profiles of cfDNA molecules of both the plasma and PMBC samples of our dataset using cfDNApro.  It requires the insert size metrics text files from Picard Collectinsertsize, and each sample should be stored in individual folder within the data_path.  The name of the folder would be the name to the sample included in the analysis.  Using the tibble data structure available by the tools, we also analysis the sample stastics and present them as a floating layer ontop of the output.

# For installation:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = '3.17')
#BiocManager::install("cfDNAPro")

#################################################################################
require(cfDNAPro)
library(scales)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(gtable)

# Set working directory
setwd("/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis")

# Get the path and list files within subdirectories
data_path <- "/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis"
list.files(data_path, full.names = TRUE,recursive = TRUE)

# Define a list for the groups/cohorts.
grp_list<-list("CAB"="CAB",
               "ROK"="ROK",
               "KSA"="KSA",
               "ACT_27"="ACT_27",
               "ROMI_31"="ROMI_31")

# Generating the plots and store them in a list.
result <-sapply(grp_list, function(x){
  result <-callSize(path = data_path) %>% 
    dplyr::filter(group==as.character(x)) %>% 
    plotSingleGroup()
}, simplify = FALSE)

# Multiplexing the plots in one figure
suppressWarnings(
  multiplex <-
    ggarrange(result$CAB$prop_plot + 
                ylim(c(0, 0.04)) +
                theme(axis.title.x = element_blank()),
              result$ROK$prop_plot + 
                theme(axis.title.x = element_blank()),
              result$KSA$prop_plot + 
                theme(axis.title.x = element_blank()),
              result$ACT_27$prop_plot + 
                theme(axis.title.x = element_blank()),
              result$ROMI_31$prop_plot + 
                theme(axis.title.x = element_blank()),
              result$CAB$cdf_plot +
                theme(axis.title.y = element_blank()),
              result$ROK$cdf_plot +
                theme(axis.title.y = element_blank()),
              result$KSA$cdf_plot +
                theme(axis.title.y = element_blank()),
              result$ACT_27$cdf_plot +
                theme(axis.title.y = element_blank()),
              result$ROMI_31$cdf_plot +
                theme(axis.title.y = element_blank()),
              labels = c("CAB", "ROK", "KSA", "ACT_27", "ROMI_31"),
              label.x = 0.2,
              ncol = 5,
              nrow = 2))

multiplex
ggsave("multiplex.png",plot = last_plot())

#  Median Size Metrics - calculate the median size fragment distribution
# Set an order for those groups (i.e. the levels of factors).
order <- c("CAB", "ROK", "KSA", "ACT_27", "ROMI_31", "healthy_control")

# Generate plots.
compare_grps<-callMetrics(data_path) %>% plotMetrics(order=order, vline=c(81,167))

# Modify plots.
p1<-compare_grps$median_prop_plot +
  ylim(c(0, 0.04)) +
  xlim(c(0, 300)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank()) +
  ggtitle("cfDNA Fragment Length in PTCL") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# based on the cumulative frequency for each insert size in a group, find the first of each which is greater than the required quartile
quartile_table <- compare_grps$median_prop_plot$data %>% group_by(group) %>%
  #filter(insert_size <= 280)  %>%
  summarize(Q1=insert_size[cdf_median >= 0.25][1],
            median=insert_size[cdf_median >= 0.5][1],
            Q3=insert_size[cdf_median >= 0.75][1],
            mean= sum(insert_size * prop_median)) %>% 
    dplyr::mutate_if(is.numeric, round, 1)

quartile_table

#add title and paddings
q_table <- tableGrob(quartile_table, rows=NULL, theme = mytheme)
title <- textGrob("Summary Statistics",gp=gpar(fontsize=12))
padding <- unit(5,"mm")
table <- gtable_add_rows(
  q_table, 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 1, 1, ncol(table))

# format text size
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.8)),
  colhead = list(fg_params=list(cex = 0.9)),
  rowhead = list(fg_params=list(cex = 0.9)))

p2<-compare_grps$median_cdf_plot +
  xlim(c(0, 300)) +
  annotation_custom(tableGrob(quartile_table, rows=NULL, theme = mytheme), 
                    xmin=20, xmax=80, ymin=0.6, ymax=0.75) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.4),
        legend.text = element_text(size = 11),
        legend.title = element_blank())

# Finalize plots.
suppressWarnings(
  median_grps<-ggpubr::ggarrange(p1,
                                 p2,
                                 label.x = 0.3,
                                 ncol = 1,
                                 nrow = 2
  ))
median_grps
ggsave("median_grps_indi.png",plot = last_plot())


## Modal Fragment Size
# Generate mode bin chart.
mode_bin <- callMode(data_path) %>% plotMode(order=order,hline = c(167,111,81))

# Show the plot.
suppressWarnings(print(mode_bin))
ggsave("mode_bin.png",plot = last_plot())

mode_stacked <- 
  callMode(data_path) %>% 
  plotModeSummary(order=order,
                  mode_partition = list(c(166,167)))
#> setting default input_type to picard.

# Modify the plot using ggplot syntax.
mode_stacked <- mode_stacked + theme(legend.position = "top")

# Show the plot.
suppressWarnings(print(mode_stacked))

## Inter-peak/valley Distance
# Plot and modify inter-peak distances.
inter_peak_dist<-callPeakDistance(path = data_path,  limit = c(50, 135)) %>%
  plotPeakDistance(order = order) +
  labs(y="Fraction") +
  theme(axis.title =  element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.5),
        legend.text = element_text(size = 11))
# Show the plot.
suppressWarnings(print(inter_peak_dist))
ggsave("inter_peak_dist.png",plot = last_plot())

# Plot and modify inter-peak distances.
inter_valley_dist<-callValleyDistance(path = data_path,  
                                      limit = c(50, 135)) %>%
  plotValleyDistance(order = order) +
  labs(y="Fraction") +
  theme(axis.title =  element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.5),
        legend.text = element_text(size = 11))
#> setting the mincount to 0. 
#>  setting the xlim to c(7,13). 
#>  setting default outfmt to df.
#> setting the mincount to 0.
#> setting default input_type to picard.

# Show the plot.
suppressWarnings(print(inter_valley_dist))
ggsave("inter_valley_dist.png",plot = last_plot())


##highlightthe peaks and valleys in the fragmentation patterns
# Calculate peaks and valleys.
peaks <- callPeakDistance(path = data_path) 
valleys <- callValleyDistance(path = data_path) 
# A line plot showing the fragmentation pattern of the example sample.
exam_plot_all <- callSize(path=data_path) %>% plotSingleGroup(vline = NULL)
# Label peaks and valleys with dashed and solid lines.
exam_plot_prop <- exam_plot_all$prop + 
  coord_cartesian(xlim = c(90,135),ylim = c(0,0.0065)) +
  geom_vline(xintercept=peaks$insert_size, colour="red",linetype="dashed") +
  geom_vline(xintercept = valleys$insert_size,colour="blue")

# Show the plot.
suppressWarnings(print(exam_plot_prop))
ggsave("exam_plot_prop.png",plot = last_plot())

# Label peaks and valleys with dots.
exam_plot_prop_dot<- exam_plot_all$prop + 
  coord_cartesian(xlim = c(90,135),ylim = c(0,0.0065)) +
  geom_point(data= peaks, 
             mapping = aes(x= insert_size, y= prop),
             color="blue",alpha=0.5,size=3) +
  geom_point(data= valleys, 
             mapping = aes(x= insert_size, y= prop),
             color="red",alpha=0.5,size=3) 
# Show the plot.
suppressWarnings(print(exam_plot_prop_dot))
ggsave("exam_plot_prop_dot.png",plot = last_plot())


##################Rerun analysis with PBMC set###########
# Define a list 
order <- c("PBMC_CAB", "PBMC_ROK", "PBMC_BHD", "PBMC_SJZ", "PBMC_LR185")

compare_grps<-callMetrics(data_path) %>% plotMetrics(order=order, vline=c(81,167))

# Modify plots.
p1<-compare_grps$median_prop_plot +
  ylim(c(0, 0.04)) +
  xlim(c(0, 300)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank()) +
  ggtitle("cfDNA Fragment Length in PTCL") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# based on the cumulative frequency for each insert size in a group, find the first of each which is greater than the required quartile
quartile_table <- compare_grps$median_prop_plot$data %>% group_by(group) %>%
  #filter(insert_size <= 280)  %>%
  summarize(Q1=insert_size[cdf_median >= 0.25][1],
            median=insert_size[cdf_median >= 0.5][1],
            Q3=insert_size[cdf_median >= 0.75][1],
            mean= sum(insert_size * prop_median)) %>% 
    dplyr::mutate_if(is.numeric, round, 1)

quartile_table

#add title and paddings
q_table <- tableGrob(quartile_table, rows=NULL, theme = mytheme)
title <- textGrob("Summary Statistics",gp=gpar(fontsize=12))
padding <- unit(5,"mm")
table <- gtable_add_rows(
  q_table, 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 1, 1, ncol(table))

# format text size
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.8)),
  colhead = list(fg_params=list(cex = 0.9)),
  rowhead = list(fg_params=list(cex = 0.9)))

p2<-compare_grps$median_cdf_plot +
  xlim(c(0, 300)) +
  annotation_custom(tableGrob(quartile_table, rows=NULL, theme = mytheme), 
                    xmin=20, xmax=80, ymin=0.6, ymax=0.75) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.4),
        legend.text = element_text(size = 11),
        legend.title = element_blank())

# Finalize plots.
suppressWarnings(
  median_grps<-ggpubr::ggarrange(p1,
                                 p2,
                                 label.x = 0.3,
                                 ncol = 1,
                                 nrow = 2
  ))
median_grps

ggsave("median_grps_PBMC.png",plot = last_plot())

