####### Script Information ########################
# Brandon P.M. Edwards
# 2021 Arboretum WBC Paper
# 1-analysis.R
# Created February 2021
# Last Updated February 2022

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
species_list <- species_list[-which(grepl("sp.", species_list, fixed = FALSE) |
                                grepl("unid.", species_list, fixed = FALSE))]

combined_red <- combined[which(combined$Species %in% species_list), ]

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
  ylim(0, max(yearly$TotalSp) + 10) +
  annotate("segment", x = 1985, xend = 1980, y = 5, yend = 11, colour = "blue") +
  annotate("segment", x = 1985, xend = 1982, y = 5, yend = 11, colour = "blue") +
  annotate("text", x = 1990, y = 3, label = "Survey Minimum: 11\n(1980 & 1982)") +
  annotate("segment", x = 2010, xend = 2021, y = 35, yend = 29, colour = "blue") +
  annotate("text", x = 2008, y = 35, label = "Survey Maximum: 29\n(2021)") +
  NULL

# Annual Shannon index
shannon_plot <- ggplot() +
  geom_line(data = yearly, aes(x = Year, y = Shannon), size = 1.25) +
  ylim(0, 3) +
  NULL

# plot species accumulation
accum_plot <- ggplot() +
  geom_line(data = yearly, aes(x = Year, y = Accumulation), size = 1.25) +
  ylim(0, max(yearly$Accumulation) + 15) +
  xlim(1980, 2029)+
  annotate("segment", x = 2024, xend = 2022, y = 45, yend = 60, colour = "blue") +
  annotate("text", x = 2024, y = 40, label = "WCSP\n(2022)") +
  annotate("segment", x = 2021, xend = 2019, y = 65, yend = 59, colour = "blue") +
  annotate("text", x = 2021, y = 70, label = "CORA*\n(2019)") +
  annotate("segment", x = 2024-10, xend = 2017, y = 40, yend = 58, colour = "blue") +
  annotate("text", x = 2014, y = 35, label = "EATO\n(2017)") +
  annotate("segment", x = 2010, xend = 2015, y = 65, yend = 57, colour = "blue") +
  annotate("text", x = 2010, y = 70, label = "WIWR\n(2015)") +
  annotate("segment", x = 2000, xend = 2013, y = 60, yend = 56, colour = "blue") +
  annotate("text", x = 1990, y = 60, label = "RBWO*, BRTH, BADO\n(2013)") +
  NULL

####### Trends #####################################

####### Select Species #############################

sp <- c("Ruffed Grouse", "Red-bellied Woodpecker", "Mourning Dove", "Evening Grosbeak")

sp_select <- combined_red[which(combined_red$Species %in% sp), ]


rugr <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Ruffed Grouse", ],
                             aes(x = Year, y = Count), alpha = 0.5) +
  geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Ruffed Grouse"),
                                      c("Year", "Count")], k = 5)),
            aes(x = Year, y = Count), size = 1.25)
rbwo <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Red-bellied Woodpecker", ],
            aes(x = Year, y = Count), alpha = 0.5) +
  geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Red-bellied Woodpecker"),
                                                 c("Year", "Count")], k = 5)),
            aes(x = Year, y = Count), size = 1.25)
modo <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Mourning Dove", ],
            aes(x = Year, y = Count), alpha = 0.5) +
  geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Mourning Dove"),
                                                 c("Year", "Count")], k = 5)),
            aes(x = Year, y = Count), size = 1.25)
evgr <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Evening Grosbeak", ],
            aes(x = Year, y = Count), alpha = 0.5) +
  geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Evening Grosbeak"),
                                                 c("Year", "Count")], k = 5)),
            aes(x = Year, y = Count), size = 1.25)

# Create matrix plot
sp_select_abundance <- ggarrange(rugr, rbwo, modo, evgr,
  nrow = 2, ncol = 2,
  labels = c("RUGR", "RBWO", "MODO", "EVGR")
)

####### Output Data and Plots #####################

png(filename = "plots/abundance.png",
    width = 6, height = 4, units = "in", res = 300)
print(yearly_abundance_plot)
dev.off()

png(filename = "plots/yearly_sp.png",
    width = 6, height = 4, units = "in", res = 300)
print(yearly_sp_plot)
dev.off()

png(filename = "plots/shannon.png",
    width = 6, height = 4, units = "in", res = 300)
print(shannon_plot)
dev.off()

png(filename = "plots/accumulation.png",
    width = 6, height = 4, units = "in", res = 300)
print(accum_plot)
dev.off()

png(filename = "plots/sp_select.png",
    width = 6, height = 6, units = "in", res = 300)
print(sp_select_abundance)
dev.off()

write.csv(combined, file = "data/combined_red.csv", row.names = FALSE)
write.csv(yearly, file = "data/yearly_metrics.csv", row.names = FALSE)
