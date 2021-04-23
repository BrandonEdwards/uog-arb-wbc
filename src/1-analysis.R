####### Script Information ########################
# Brandon P.M. Edwards
# 2021 Arboretum WBC Paper
# 1-analysis.R
# Created February 2021
# Last Updated February 2021

####### Import Libraries and External Files #######

library(zoo)
library(ggplot2)
library(GGally)
library(ggpubr)
theme_set(theme_pubclean())

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

####### Yearly Abundance, Richness, Diversity #####

# Get rid of non-species
combined_red <- combined[which(combined$Species %in% summary_stats[["Species"]]), ]

# Total abundance by year and rolling mean
yearly <- aggregate(Count ~ Year, data = combined_red, FUN = sum)
names(yearly) <- c("Year", "Count")
yearly_abundance_rolling <- data.frame(rollmean(yearly, k = 5))
names(yearly_abundance_rolling) <- c("Year", "Rolling")
yearly <- merge(yearly, yearly_abundance_rolling, by = "Year", all = TRUE)

# Total species each year
yearly <- merge(yearly, count_data[, c("Year", "TotalSp")], by = "Year")

# Calculate Shannon index for each year
yearly$Shannon <- NA
for (y in unique(yearly$Year))
{
  temp <- combined_red[which(combined_red$Year == y &
                               combined_red$Count > 0),]
  p_vec <- temp$Count/yearly[yearly$Year == y, "Count"]
  yearly[yearly$Year == y, "Shannon"] <- -1 * (sum(p_vec * log(p_vec)))
}

# Calculate species accumulation each year
yearly$Accumulation <- NA
species_discovered <- NULL
for (y in unique(yearly$Year))
{
  temp <- combined_red[which(combined_red$Year == y &
                               combined_red$Count > 0),]
  if (length(species_discovered) == 0)
  {
    species_discovered <- temp$Species
  }else
  {
    species_discovered <- union(species_discovered, temp$Species)
  }
  
  yearly[yearly$Year == y, "Accumulation"] <- length(species_discovered)
}

# Plot yearly abundance and rolling mean
yearly_abundance_plot <- ggplot() +
  geom_line(data = yearly, aes(x = Year, y = Count), alpha = 0.5) +
  geom_line(data = yearly, aes(x = Year, y = Rolling), size = 1.25) +
  NULL

# plot total species by year
yearly_sp_plot <- ggplot() +
  geom_line(data = yearly, aes(x = Year, y = TotalSp), size = 1.25) +
  ylim(0, max(yearly$TotalSp) + 5) +
  NULL

# Annual Shannon index
shannon_plot <- ggplot() +
  geom_line(data = yearly, aes(x = Year, y = Shannon), size = 1.25) +
  ylim(0, 3) +
  NULL

# plot species accumulation
accum_plot <- ggplot() +
  geom_line(data = yearly, aes(x = Year, y = Accumulation), size = 1.25) +
  ylim(0, max(yearly$Accumulation) + 5) +
  NULL
####### Output Data ###############################

#write.csv(combined, file = "data/combined.csv", row.names = FALSE)
