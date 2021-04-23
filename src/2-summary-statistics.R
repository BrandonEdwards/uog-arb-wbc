####### Script Information ########################
# Brandon P.M. Edwards
# 2021 Arboretum WBC Paper
# 2-summary-statistics.R
# Created February 2021
# Last Updated February 2021

####### Import Libraries and External Files #######


####### Read Data #################################

arb <- read.csv("data/combined.csv")

####### Species Diversity #########################

summary_stats <- vector(mode = "list", length = 0)

sp_list <- unique(arb$Species)
sp_list_rem <- sp_list[-which(grepl("sp.", sp_list, fixed = FALSE) |
                                grepl("unid.", sp_list, fixed = FALSE))]

summary_stats[["Species"]] <- sp_list_rem
summary_stats[["Total_Species"]] <- length(sp_list_rem)

####### Output Data ###############################

save(summary_stats, file = "data/summary_stats.rda")
