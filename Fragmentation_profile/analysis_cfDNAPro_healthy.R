# Authors: Friday
# Modified Date: 30-05-2023
# Version: 1.0

#################################################################################
# This R script is designed to produce fragmentation profiles for healthy controls
# and analyse the distribution stastics.
#################################################################################
library(cfDNAPro)
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
data_path <- "/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis/healthy"
list.files(data_path, full.names = TRUE,recursive = TRUE)

#  Median Size Metrics - calculate the median size fragment distribution
# Set an order for those groups (i.e. the levels of factors).
order <- c("Ctrl_1","Ctrl_2","Ctrl_3","Ctrl_4","Ctrl_5","Ctrl_6","Ctrl_7",
           "Ctrl_8","Ctrl_9","Ctrl_10","Ctrl_11","Ctrl_12","Ctrl_13","Ctrl_14",
           "Ctrl_15","Ctrl_16","Ctrl_17","Ctrl_18","Ctrl_19","Ctrl_20","Ctrl_21",
           "Ctrl_22")

# Generate plots.
compare_grps<-callMetrics(data_path) %>% plotMetrics(order=order, vline=c(81,167))

# based on the cumulative frequency for each insert size in a group, find the first of each which is greater than the required quartile
quartile_table <- compare_grps$median_prop_plot$data %>% group_by(group) %>%
  summarize(Q1=insert_size[cdf_median >= 0.25][1],
            median=insert_size[cdf_median >= 0.5][1],
            Q3=insert_size[cdf_median >= 0.75][1],
            mean= sum(insert_size * prop_median)) %>% 
    dplyr::mutate_if(is.numeric, round, 1)

quartile_table

# Modify plots.
p1<-compare_grps$median_prop_plot +
  ylim(c(0, 0.04)) +
  xlim(c(0, 500)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold")) +
  theme(legend.position = c(0.85, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank()) +
  ggtitle("cfDNA Fragment Length in healthy controls") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2<-compare_grps$median_cdf_plot +
  xlim(c(0, 500)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.position = c(0.85, 0.5),
        legend.text = element_text(size = 11),
        legend.title = element_blank())

suppressWarnings(
  median_grps<-ggpubr::ggarrange(p1,
                                 p2,
                                 label.x = 0.3,
                                 ncol = 1,
                                 nrow = 2
  ))
median_grps
ggsave("median_grps_healthy.png",plot = last_plot())


##Compare Group 20 and 19 which has the largest difference
order <- c("Ctrl_20","Ctrl_19")
compare_grps<-callMetrics(data_path) %>% plotMetrics(order=order, vline=c(81,167))
p1<-compare_grps$median_prop_plot +
  ylim(c(0, 0.04)) +
  xlim(c(0, 500)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold")) +
  theme(legend.position = c(0.85, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank()) +
  ggtitle("cfDNA Fragmentations in healthy controls_20 & 19") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2<-compare_grps$median_cdf_plot +
  xlim(c(0, 500)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.position = c(0.85, 0.5),
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
ggsave("median_grps_healthy_max_dif.png",plot = last_plot())

#perform Two-sample Kolmogorov-Smirnov test
data1 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(cur_group()$group == "Ctrl_20")
data2 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(cur_group()$group == "Ctrl_19")
ks.test(data1$cdf_median, data2$cdf_median)


#perform Two-sample Kolmogorov-Smirnov test with cut off 200 insert size
data1 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(insert_size <= 200)  %>%  filter(cur_group()$group == "Ctrl_20")
data2 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(insert_size <= 200)  %>%  filter(cur_group()$group == "Ctrl_19")
ks.test(data1$cdf_median, data2$cdf_median)
# D = 0.1345, p-value = 0.09067

#perform Two-sample Kolmogorov-Smirnov test with cut off 180 insert size
data1 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(insert_size <= 180)  %>%  filter(cur_group()$group == "Ctrl_20")
data2 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(insert_size <= 180)  %>%  filter(cur_group()$group == "Ctrl_19")
ks.test(data1$cdf_median, data2$cdf_median)
# D = 0.072848, p-value = 0.8178

#perform Two-sample Kolmogorov-Smirnov test with cut off 81-200 insert size
data1 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(insert_size <= 200)  %>%  filter(insert_size >= 81)  %>%  filter(cur_group()$group == "Ctrl_20")
data2 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(insert_size <= 200)  %>%  filter(insert_size >= 81) %>%  filter(cur_group()$group == "Ctrl_19")
ks.test(data1$cdf_median, data2$cdf_median)
#D = 0.19167, p-value = 0.02435

