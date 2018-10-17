# This file compares the two groups directly with each other


# 1) Preliminaries

# define file path
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# load packages
library(TraMineR)
library(viridis)
library(cluster)
library(WeightedCluster)

# load transformed data
matches <- read.csv(paste0(path, "matches.csv"))
treat   <- subset(matches, matches$group == T)
control <- subset(matches, matches$group == F)

rm(matches)

# 2) defining the sequence object

# setting sequence labels 
longlab <- c("Single", "Single, Children", 
             "Cohabitation", "Cohabitation, Children", 
             "Marriage", "Marriage, Children")
shortlab <- c("SNC", "SC", "CNC", "CC", "MNC", "MC")

# defining custom legend colors
custom_col <- viridis(n = length(shortlab), begin = 1, end = 0)

# selecting sequence columns
seqcols <- grep("state.", names(treat))

# defining the sequences
treat_seq   <- seqdef(treat, var = seqcols, states = shortlab, labels = longlab,
                 cpal = custom_col)
control_seq <- seqdef(control, var = seqcols, states = shortlab,
                    labels = longlab, cpal = custom_col)

# defining the cost matrices
treat_dist   <- seqdist(treat_seq, method = "OM", indel = 2, sm = "CONSTANT")
control_dist <- seqdist(control_seq, method = "OM", indel = 2, sm = "CONSTANT")

# 3) Clustering

# 3.1) first, for treatment
treat_clust <- hclust(as.dist(treat_dist), method = "ward.D")

# evaluation of cluster sizes: 5 and 4 (as before)
summary(as.clustrange(treat_clust, diss = treat_dist, ncluster = 10), max.rank = 4)

# 5 clusters:
om5_treat <- cutree(treat_clust, k = 5)
om5_treat <- factor(om5_treat, levels = c(2, 3, 1, 4, 5),
                  labels = c("Early Family Formation", "Later Family Formation", 
                             "Childless Marriage", "Childless Individuals", 
                             "Non-traditional Parenthood"))
# seqdplot(treat_seq, group = om5_treat, xt = 15:60)

# 4 clusters:
om4_treat <- cutree(treat_clust, k = 4)
om4_treat <- factor(om4_treat, levels = c(2, 3, 1, 4),
                  labels = c("Early Family Formation", "Later Family Formation", 
                             "Childless Marriage", "Non-traditional Family Form"))
# seqdplot(treat_seq, group = om4_treat, xt = 15:60)

# 3.2) Now, cluster control group
control_clust <- hclust(as.dist(control_dist), method = "ward.D")

# evaluation of cluster sizes: 5 and 4 (same as for treatment)
summary(as.clustrange(control_clust, diss = control_dist, ncluster = 10), max.rank = 4)

# 5 clusters:
om5_control <- cutree(control_clust, k = 5)
om5_control <- factor(om5_control, levels = c(5, 3, 2, 1, 4),
                    labels = c("Early Family Formation", "Later Family Formation", 
                               "Childless Marriage", "Childless Individuals", 
                               "Non-traditional Parenthood"))
seqfplot(control_seq, group = om5_control, xt = 15:60)

# 4 clusters:
om4_control <- cutree(control_clust, k = 4)
om4_control <- factor(om4_control, levels = c(4, 3, 2, 1),
                    labels = c("Early Family Formation", "Later Family Formation", 
                               "Childless Marriage", "Non-traditional Family Form"))
seqdplot(control_seq, group = om4_control, xt = 15:60)

# for pairswise 1:1, same four 
# for propensity 1:1, same four
# for mahalanobis 1:1, same four