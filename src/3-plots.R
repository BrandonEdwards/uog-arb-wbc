####### Script Information ########################
# Brandon P.M. Edwards
# 2021 Arboretum WBC Paper
# 3-plots.R
# Created February 2021
# Last Updated February 2021

####### Import Libraries and External Files #######

library(zoo)
library(ggplot2)
library(GGally)
library(ggpubr)
theme_set(theme_pubclean())

####### Read Data #################################

arb <- read.csv("data/combined.csv")
sample <- read.csv("data/count_data.csv")
load("data/summary_stats.rda")

####### Total Abundance and Species by Year########

# Get rid of non-species
arb_red <- arb[which(arb$Species %in% summary_stats[["Species"]]), ]

# Total abundance by year
yearly_abundance <- aggregate(Count ~ Year, data = arb_red, FUN = sum)
names(yearly_abundance) <- c("Year", "Count")

yearly_abundance_rolling <- rollmean(yearly_abundance, k = 5)

yearly_abundance_plot <- ggplot() +
  geom_line(data = yearly_abundance, aes(x = Year, y = Count), alpha = 0.5) +
  geom_line(data = data.frame(yearly_abundance_rolling), aes(x = Year, y = Count), size = 1.25) +
  #geom_line(data = arb_red, aes(x = Year, y = Count, group = Species),alpha = 0.2) +
  NULL

# Total species by year
yearly_sp_plot <- ggplot() +
  geom_line(data = sample, aes(x = Year, y = TotalSp))

####### Output Data ###############################

write.csv(combined, file = "data/combined.csv", row.names = FALSE)
