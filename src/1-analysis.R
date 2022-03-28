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
library(lubridate)
library(stringr)
library(magrittr)
library(reshape)
theme_set(theme_pubclean())

####### Read Data #################################

count_data <- read.csv("data/count_data.csv")
bird_data <- read.csv("data/bird_data.csv")
cbc_data <- read.csv("data/ONGU-transformed-per-hour.csv")

####### Wrangle Data ##############################

# Melt CBC data
cbc_data <- melt(cbc_data, id = "Species")
names(cbc_data) <- c("Species", "Year", "BirdsPerHour")
cbc_data$Year <- as.numeric(substr(cbc_data$Year, start = 2, stop = 6))

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

# Get rid of non-species
species_list <- species_list[-which(grepl("sp.", species_list, fixed = FALSE) |
                                      grepl("unid.", species_list, fixed = FALSE))]

combined_red <- combined[which(combined$Species %in% species_list), ]

####### Calculate Survey Length & Birds Per Hour ##

combined_red$TimeStart <- str_pad(combined_red$TimeStart,
                                  4,
                                  pad = "0")
combined_red$TimeStart <- format(strptime(sprintf('%s', combined_red$TimeStart), format='%H%M'), '%H:%M')
  

combined_red$TimeEnd <- str_pad(combined_red$TimeEnd,
                                  4,
                                  pad = "0")
combined_red$TimeEnd <- format(strptime(sprintf('%s', combined_red$TimeEnd), format='%H%M'), '%H:%M')

combined_red$SurveyLength <- as.numeric(hm(combined_red$TimeEnd) - hm(combined_red$TimeStart)) / 3600

combined_red$BirdsPerHour <- combined_red$Count / combined_red$SurveyLength

####### Yearly Abundance, Richness, Diversity #####

# Total abundance by year and rolling mean
yearly <- aggregate(BirdsPerHour ~ Year, data = combined_red, FUN = sum)
names(yearly) <- c("Year", "BirdsPerHour")
yearly_abundance_rolling <- data.frame(rollmean(yearly, k = 5))
names(yearly_abundance_rolling) <- c("Year", "Rolling")
yearly <- merge(yearly, yearly_abundance_rolling, by = "Year", all = TRUE)

# Total species each year
yearly <- merge(yearly, count_data[, c("Year", "TotalSp")], by = "Year")

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
  ylim(0, 800)+
  geom_line(data = yearly, aes(x = Year, y = Count), size = 1.25) +#alpha = 0.5) +
  #geom_line(data = yearly, aes(x = Year, y = Rolling), size = 1.25) +
  annotate("segment", x = 1996, xend = 1987, y = 42, yend = 42, colour = "blue") +
  annotate("text", x = 2004, y = 42, label = "Survey Minimum: 42 (1987)") +
  annotate("segment", x = 1996+6, xend = 1987+6, y = 733, yend = 733, colour = "blue") +
  annotate("text", x = 2004+6, y = 733, label = "Survey Maximum: 733\n(1993)") +
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
  annotate("text", x = 1990, y = 60, label = "RBWO*, BRTH*, BADO\n(2013)") +
  NULL

####### Gulls ######################################

gulls <- c("Herring Gull", "Ring-billed Gull", "Great Black-backed Gull",
           "Glaucous Gull")

gull_select <- combined_red[which(combined_red$Species %in% gulls), ]
gull_select$Species <- factor(gull_select$Species,
                              levels = c("Ring-billed Gull",
                                         "Herring Gull",
                                         "Great Black-backed Gull",
                                         "Glaucous Gull"))

gull_plot <- ggplot(data = gull_select) +
  geom_line(aes(x = Year, y = BirdsPerHour, group = Species, color = Species)) +
#  stat_summary(aes(x = Year, y = Count), fun = "sum", geom = "line", size = 1.25) +
  theme(legend.position = "right") +
  #ylim(0, 450) +
  annotate("segment", x = 2003, xend = 2003, y = 0, yend = 100, colour = "blue", size = 1.5) +
  annotate("text", x = 2012, y = 50, label = "Guelph Landfill\nClosure (2003)") +
  NULL

####### Select Species #############################

sp <- c("Ruffed Grouse", "Red-bellied Woodpecker", "Mourning Dove", "Evening Grosbeak")

sp_select <- combined_red[which(combined_red$Species %in% sp), ]; sp_select$dataset <- "WBC"
cbc_select <- cbc_data[which(cbc_data$Species %in% sp), ]; cbc_select$dataset <- "CBC"

wbc_cbc <- rbind(sp_select[, c("Species", "Year", "BirdsPerHour", "dataset")],
                 cbc_select[, c("Species", "Year", "BirdsPerHour", "dataset")])
wbc_cbc <- wbc_cbc[which(wbc_cbc$Year >= 1980 & wbc_cbc$Year <= 2020), ]
wbc_cbc <- wbc_cbc[-which(wbc_cbc$Year == 1983), ]

corr_coefs <- vector(mode = "list", length = length(sp))
names(corr_coefs) <- sp
for (s in sp)
{
  temp <- wbc_cbc[which(wbc_cbc$Species == s), ]
  corr_coefs[[s]] <- cor(temp[temp$dataset == "WBC", "BirdsPerHour"],
                         temp[temp$dataset == "CBC", "BirdsPerHour"])
}


rugr <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Ruffed Grouse", ],
                             aes(x = Year, y = BirdsPerHour), size = 1.25) +
  geom_line(data = cbc_select[cbc_select$Species == "Ruffed Grouse", ], 
            aes(x = Year, y = BirdsPerHour)) +
  annotate("text", x = 2010, y = 2.5, label = paste0("r = ", round(corr_coefs[["Ruffed Grouse"]], digits = 2))) +
 # geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Ruffed Grouse"),
  #                                    c("Year", "BirdsPerHour")], k = 5)),
           # aes(x = Year, y = BirdsPerHour), size = 1.25) +
  ylim(0, 3) +
  xlim(1980, 2021) +
  ylab("Birds Per Party Hour") +
  NULL

rbwo <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Red-bellied Woodpecker", ],
            aes(x = Year, y = BirdsPerHour), size = 1.25) +
  geom_line(data = cbc_select[cbc_select$Species == "Red-bellied Woodpecker", ], 
            aes(x = Year, y = BirdsPerHour)) +
  annotate("text", x = 2010, y = 2.5, label = paste0("r = ", round(corr_coefs[["Red-bellied Woodpecker"]], digits = 2))) +
  #geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Red-bellied Woodpecker"),
   #                                              c("Year", "BirdsPerHour")], k = 5)),
            #aes(x = Year, y = BirdsPerHour), size = 1.25) +
  ylim(0,3) +
  xlim(1980, 2021) +
  ylab("Birds Per Party Hour") +
  NULL

modo <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Mourning Dove", ],
            aes(x = Year, y = BirdsPerHour), size = 1.25) +
  geom_line(data = cbc_select[cbc_select$Species == "Mourning Dove", ], 
            aes(x = Year, y = BirdsPerHour)) +
  annotate("text", x = 2010, y = 17, label = paste0("r = ", round(corr_coefs[["Mourning Dove"]], digits = 2))) +
  #geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Mourning Dove"),
   #                                              c("Year", "BirdsPerHour")], k = 5)),
            #aes(x = Year, y = BirdsPerHour), size = 1.25) +
  ylim(0,20) +
  xlim(1980, 2021) +
  ylab("Birds Per Party Hour") +
  NULL

evgr <- ggplot() + 
  geom_line(data = sp_select[sp_select$Species == "Evening Grosbeak", ],
            aes(x = Year, y = BirdsPerHour), size = 1.25) +
  geom_line(data = cbc_select[cbc_select$Species == "Evening Grosbeak", ], 
            aes(x = Year, y = BirdsPerHour)) +
  annotate("text", x = 2010, y = 17, label = paste0("r = ", round(corr_coefs[["Evening Grosbeak"]], digits = 2))) +
  #geom_line(data = data.frame(rollmean(sp_select[which(sp_select$Species == "Evening Grosbeak"),
   #                                              c("Year", "BirdsPerHour")], k = 5)),
           # aes(x = Year, y = BirdsPerHour), size = 1.25) +
  ylim(0,20) +
  xlim(1980, 2021) +
  ylab("Birds Per Party Hour") +
  NULL

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

png(filename = "plots/gulls.png",
    width = 6, height = 4, units = "in", res = 300)
print(gull_plot)
dev.off()

write.csv(combined, file = "data/combined_red.csv", row.names = FALSE)
write.csv(yearly, file = "data/yearly_metrics.csv", row.names = FALSE)
