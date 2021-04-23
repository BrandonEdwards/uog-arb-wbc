####### Script Information ########################
# Brandon P.M. Edwards
# 2021 Arboretum WBC Paper
# 1-combine-data.R
# Created February 2021
# Last Updated February 2021

####### Import Libraries and External Files #######


####### Read Data #################################

count_data <- read.csv("data/count_data.csv")
bird_data <- read.csv("data/bird_data.csv")

####### Wrangle Data ##############################

species_list <- unique(bird_data$Species)

combined <- vector(mode = "list", length = length(species_list))
names(combined) <- species_list

for(sp in species_list)
{
  bird_sp <- bird_data[which(bird_data$Species == sp), ]
  combined_temp <- merge(bird_sp, count_data, by = "Year", all = TRUE)
  combined_temp$Species <- sp
  combined_temp$Count <- as.numeric(as.character(combined_temp$Count))
  combined_temp$Count[is.na(combined_temp$Count)] <- 0
  
  combined[[sp]] <- combined_temp
}

combined <- do.call(rbind, combined)

drops <- c("TotalInd", "TotalSp")
combined <- combined[, !(names(combined) %in% drops)]

####### Output Data ###############################

write.csv(combined, file = "data/combined.csv", row.names = FALSE)
