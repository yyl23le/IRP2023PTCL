require(cfDNAPro)
library(scales)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(gtable)

# Set working directory
setwd("/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis/SNP/CAB_FAS")

# Get the path and list files within subdirectories
data_path <- "/home/yyl23/Desktop/BS7130_Independent_Research_Project/analysis/SNP/CAB_FAS_JV"
list.files(data_path, full.names = TRUE,recursive = TRUE)

#  Median Size Metrics - calculate the median size fragment distribution
# Set an order for those groups (i.e. the levels of factors).
order <- c("alt", "ref","unfiltered")

# Generate plots.
compare_grps<-callMetrics(data_path) %>% plotMetrics(order=order, vline=c(81,167))

# based on the cumulative frequency for each insert size in a group, find the first of each which is greater than the required quartile
quartile_table <- compare_grps$median_prop_plot$data %>% group_by(group) %>%
  #filter(insert_size <= 280)  %>%
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
  ggtitle("cfDNA Fragment Length in PTCL") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# format text size
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.8)),
  colhead = list(fg_params=list(cex = 0.9)),
  rowhead = list(fg_params=list(cex = 0.9)))

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

p2<-compare_grps$median_cdf_plot +
  xlim(c(0, 500)) +
  annotation_custom(table, 
                    xmin=50, xmax=80, ymin=0.7, ymax=0.85) +
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

ggsave("median_grps_cab_fas_with_control.png",plot = last_plot())

#perform Two-sample Kolmogorov-Smirnov test
data1 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(cur_group()$group == "alt")
data2 <- compare_grps$median_prop_plot$data %>% group_by(group) %>%  filter(cur_group()$group == "ref")
ks.test(data1$cdf_median, data2$cdf_median)
