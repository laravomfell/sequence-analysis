# This file plots and cluster the family formation sequences of WWII survivors

# It consists of XXX parts
# 1. defining the WWII survivor sequences
# 2. defining the sequence distances: equal, timing, order
# 3. optimal matching clustering and evaluating the clusters
# 4. clustering beyond optimal matching (centroid, lcs, lcp)
# 5. evaluating cluster membership

# --------------
# Preliminaries
# --------------

# set seed for reproducibility of cluster labels
set.seed(123)

# define file path
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# load packages
library(TraMineR)
library(viridis)
library(cluster)
library(WeightedCluster)
library(ggplot2)

# load transformed data
formation <- read.csv(paste0(path, "formation.csv"))

# --------------------------------
# 
# 1. defining the sequence object
#
# --------------------------------

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

# -----------------------------------
#
# 2. defining the distance measures
#
# -----------------------------------

# define three different distance measures

# 2.1 equal important of timing and order
omdist_equal <- seqdist(WW_seq, method = "OM", indel = 2, sm = "CONSTANT")

# 2.2 timing over order
omdist_timing <- seqdist(WW_seq, method = "OM", indel = 4, sm = "CONSTANT")

# 2.3 order over timing
omdist_order <- seqdist(WW_seq, method = "OM", indel = 99, sm = ccost)

# ----------------------------
#
# 3. Clustering the sequences
#
# ----------------------------

# 3.1 equal cost metric

# cluster
WW_clust_equal <- hclust(as.dist(omdist_equal), method = "ward.D")

# evaluation: testing different cluster sizes
WW_clust_equal_range <- as.clustrange(WW_clust_equal, 
                                      diss = omdist_equal, ncluster = 10)
summary(WW_clust_equal_range, max.rank = 4)
plot(WW_clust_equal_range, 
     stat = c("ASW", "CHsq", "R2"), 
     norm = "zscore")

# Interpretation: cluster sizes of 5 and 4

# ----------
# Size five
# ----------
om5 <- cutree(WW_clust_equal, k = 5)

# Interpreting the clusters
om5info <- factor(om5, levels = c(2, 3, 1, 4, 5),
                  labels = c("Early Family Formation", "Later Family Formation", 
                             "Childless Marriage", "Childless Individuals", 
                             "Non-traditional Parenthood"))
# state distribution plot
seqdplot(WW_seq, group = om5info, xt = 15:60)
# sequence frequency plot
seqfplot(WW_seq, group = om5info, xt = 15:60)
# mean time plot
seqmtplot(WW_seq, group = om5info)

# ----------
# Size four
# ----------
om4 <- cutree(WW_clust_equal, k = 4)

# Interpretation
om4info <- factor(om4, levels = c(2, 3, 1, 4),
                  labels = c("Early Family Formation", "Later Family Formation", 
                             "Childless Marriage", "Non-traditional Family Form"))
# state distribution plot
seqdplot(WW_seq, group = om4info, xt = 15:60)
# sequence frequency plot
seqfplot(WW_seq, group = om4info, xt = 15:60)
# mean time plot
seqmtplot(WW_seq, group = om4info)


# 3.2 timing over order

# cluster
WW_clust_timing <- hclust(as.dist(omdist_timing), method = "ward.D")

# evaluate
WW_clust_timing_range <- as.clustrange(WW_clust_timing, diss = omdist_timing, ncluster = 10)
plot(WW_clust_equal_range, 
     stat = c("ASW", "CHsq", "R2"), 
     norm = "zscore")
# Interpretation: same cluster sizes (5, 4)

# Preliminary result: yields the exact same clusters
# briefly illustrated here but skipped for brevity
om5time <- cutree(WW_clust_timing, k = 5)
seqdplot(WW_seq, group = om5time, xt = 15:60)


# 3.3 order over timing

# cluster
WW_clust_timing <- hclust(as.dist(omdist_timing), method = "ward.D")

# evaluate
WW_clust_timing_range <- as.clustrange(WW_clust_timing, diss = omdist_timing, ncluster = 10)
plot(WW_clust_timing_range, 
     stat = c("ASW", "CHsq", "R2"), 
     norm = "zscore")
# same as 3.2

# ----------------------------
#
# 4. other clustering methods
#
# ----------------------------


# 4.1 Centroid Clustering
clustercent <- hclust(as.dist(omdist_equal), method = "centroid")

# evaluate
centrange <- as.clustrange(clustercent, diss = omdist_equal, ncluster = 10)
summary(centrange, max.rank = 4)
plot(centrange, stat = c("ASW", "CHsq", "R2"), norm = "zscore")
# Interpretation: up to 8-10 clusters

# given the small sample size, consider eg 6 clusters
cent6 <- cutree(clustercent, k = 6)
seqdplot(WW_seq, group = cent6)
# three of the six clusters contain less than five sequences
# -> do not use centroid clustering

# 4.2 longest common prefix
lcpdist <- seqdist(WW_seq, method = "LCP")
WW_clust_lcp <- hclust(as.dist(lcpdist), method = "ward.D")

# Evaluation
WW_clust_lcp_range <- as.clustrange(WW_clust_lcp, diss = lcpdist, ncluster = 10)
summary(WW_clust_lcp_range, max.rank = 4)
plot(WW_clust_lcp_range, stat = c("ASW", "CHsq", "R2"), norm = "zscore")
# Interpretation: suggests 7-8 clusters

lcp7 <- cutree(WW_clust_lcp, k = 7)
seqdplot(WW_seq, group = lcp7) 
# very inconclusive clusters

# 4.3 longest common subsequence
lcsdist <- seqdist(WW_seq, method = "LCS")
WW_clust_lcs <- hclust(as.dist(lcsdist), method = "ward.D")


# Evaluation
WW_clust_lcs_range <- as.clustrange(WW_clust_lcs, diss = lcsdist, ncluster = 10)
plot(WW_clust_lcs_range, stat = c("ASW", "CHsq", "R2"), norm = "zscore")
# Interpretation: suggests 5 clusters

lcs5 <- cutree(WW_clust_lcs, k = 5)
seqdplot(WW_seq, group = lcs5)
# same clusters as in optimal matching

# ----------------------
#
# 4. Cluster membership
#
# ---------------------

# cluster membership per country
round(prop.table(table(formation$country, om5info), margin = 1), 2)

# compare membership between birth cohorts
formation$WW_age_cohort <- cut(formation$WW_age, c(0, 5, 10, 15, 20, 35), include.lowest = T)

round(prop.table(table(formation$WW_age_cohort, om5info), margin = 1), 2)

# plot the cohorts

# define helper dataset for plotting
cohort_helperframe <- ftable(prop.table(table(formation$WW_age_cohort, om4info), margin = 1))
cohort_helperframe <- as.data.frame(cohort_helperframe)

# create cohort labels
coh <- c("0-5", "6-10", "11-15", "16-20", "21-35")
cohort_helperframe$coh <- factor(coh, levels = c("0-5", "6-10", "11-15", "16-20", "21-35"))

# create bar chart of the proportions
ggplot(cohort_helperframe, aes(x = om4info, y = Freq, fill = coh)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_viridis(discrete = T) +
  ylab("Frequency") +
  xlab("Clusters") + 
  labs(fill = "WWII Cohort") +
  theme_bw()