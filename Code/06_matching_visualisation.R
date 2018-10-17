# define file path
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# load packages
library(TraMineR)
library(viridis)
library(TraMineRextras)

# load transformed data
matches <- read.csv(paste0(path, "matches2.csv"))

# 2) defining the sequence object

# setting sequence labels 
longlab <- c("Single", "Single, Children", 
             "Cohabitation", "Cohabitation, Children", 
             "Marriage", "Marriage, Children")
shortlab <- c("SNC", "SC", "CNC", "CC", "MNC", "MC")

# defining custom legend colors
custom_col <- viridis(n = length(shortlab), begin = 1, end = 0)

# defining the sequence
seqcols <- grep("state.", names(matches))
match_seq <- seqdef(matches, var = seqcols, states = shortlab, labels = longlab,
                    cpal = custom_col)
nonWW_seq <- seqdef(matches[matches$group == F,], var = seqcols, states = shortlab,
                    labels = longlab, cpal = custom_col)

# Plot 1: Sequence Index Plot, unsorted
seqIplot(nonWW_seq, xt = 15:60, 
         xlab = "Age", ylab = "N Sequences", 
         sortv = matches$marage[matches$group == F],
         with.legend = F)

# Plot 2: Sequence Index Plot, sorted by gender
seqIplot(nonWW_seq, xt = 15:60, 
         xlab = "Age", ylab = "N Sequences", 
         sortv = matches$marage[matches$group == F],
         group = matches$gender[matches$group == F],
         with.legend = F)

# Plot 3: Sequence Frequency Plot
seqfplot(nonWW_seq, xt = 15:60, xlab = "Age", with.legend = F)

# Plot 4: Sequence State distribution
seqdplot(nonWW_seq, xt = 15:60, ylab = "Percent of Sequences", with.legend = F)

seqplot.rf(nonWW_seq, diss = omdist_equal, k = 100, sortv = matches$marage[matches$group == F], xtlab = 15:60)
