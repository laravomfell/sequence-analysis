# This file visualizes the complete sequences only. No clustering.

# 1) Preliminaries

# define file path
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# picture path
picturepath <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Paper1"

# load necessary library for plotting
library(ggplot2)
library(prettyunits)
library(TraMineR)
library(viridis)
library(TraMineRextras)

# load transformed data
formation <- read.csv(paste0(path, "formation.csv"))

# Plot 1: distribution of age at end of WWII per gender

ggplot(formation, aes(x = WW_age)) + 
  geom_histogram(binwidth=1) + 
  facet_grid(~gender) +
  labs(x = "Age in 1944", y = "Count") +
  scale_x_continuous(breaks = pretty(formation$WW_age))
# uncomment if saving picture
# ggsave(paste0(picturepath, "/freq_age_gender.png"))

# Table 1: country per gender and total
freq <- as.data.frame(cbind(xtabs( ~ country + gender, data = formation)))
freq$Total <- freq$Female + freq$Male

# PLOTTING THE SEQUENCES


# 2) setting the sequences

# setting sequence labels 
longlab <- c("Single", "Single, Children", 
             "Cohabitation", "Cohabitation, Children", 
             "Marriage", "Marriage, Children")
shortlab <- c("SNC", "SC", "CNC", "CC", "MNC", "MC")

# defining custom legend colors
custom_col <- viridis(n = length(shortlab), begin = 1, end = 0)

# defining the sequence
seqcols <- grep("state.", names(formation))
WW_seq <- seqdef(formation, var = seqcols, states = shortlab, labels = longlab,
                 cpal = custom_col)

# 3) Plotting the sequences

# Plot 1: Sequence Index Plot, unsorted
seqIplot(WW_seq, xt = 15:60, 
         xlab = "Age", ylab = "N Sequences", 
         sortv = formation$marage,
         with.legend = F)

# Plot 2: Sequence Index Plot, sorted by gender
seqIplot(WW_seq, xt = 15:60, 
         xlab = "Age", ylab = "N Sequences", 
         sortv = formation$marage,
         group = formation$gender,
         with.legend = F)

# Plot 3: Sequence Frequency Plot
seqfplot(WW_seq, xt = 15:60, xlab = "Age", with.legend = F)

# Plot 4: Sequence State distribution
seqdplot(WW_seq, xt = 15:60, ylab = "Percent of Sequences", with.legend = F)

# Plot 5: Mean Time Plot
seqmtplot(WW_seq, ylab = "Average Number of Years")


# Relative Frequency Sequence Plot
seqplot.rf(WW_seq, diss = omdist_equal, k = 100, sortv = formation$marage, xtlab = 15:60)
# my interpretation: for many people, very stable but consistent outliers across all medoids


# some ideas about entropy
formation$entropy <- seqient(WW_seq)
formation$cohort <- cut(formation$yrbirth, c(1910, 1918, 1928, 1938, 1944), include.lowest = T)
boxplot(entropy ~ cohort, data = formation)

