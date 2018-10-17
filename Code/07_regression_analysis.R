# first create clusters for the entire group

# testing multinomial cluster membership for the groups

set.seed(123)
# 1) Preliminaries

# define file path
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# load packages
library(TraMineR)
library(viridis)
library(cluster)
library(WeightedCluster)
library(nnet)
library(plyr)
library(dplyr)

# load transformed data, preserve matching classes
matches <- read.csv(paste0(path, "matches.csv"), colClasses = c(rep(NA, 57), "character"))

# 2) defining the sequence object

# setting sequence labels 
longlab <- c("Single", "Single, Children", 
             "Cohabitation", "Cohabitation, Children", 
             "Marriage", "Marriage, Children")
shortlab <- c("SNC", "SC", "CNC", "CC", "MNC", "MC")

# defining custom legend colors
custom_col <- viridis(n = length(shortlab), begin = 1, end = 0)

# selecting sequence columns
seqcols <- grep("state.", names(matches))

# defining the sequences
match_seq <- seqdef(matches, var = seqcols, states = shortlab, 
                    labels = longlab, cpal = custom_col)

# defining the cost matrices
match_dist <- seqdist(match_seq, method = "OM", indel = 2, sm = "CONSTANT")

# 3) Clustering

match_clust <- hclust(as.dist(match_dist), method = "ward.D")

# evaluation of cluster sizes: 5 and 4 (as before)
summary(as.clustrange(match_clust, diss = match_dist, ncluster = 10), max.rank = 4)

# 5 clusters:
om5_match <- cutree(match_clust, k = 5)
om5_match <- factor(om5_match, levels = c(4, 3, 2, 1, 5),
                    labels = c("Early Family Formation", "Later Family Formation", 
                               "Childless Marriage", "Childless Individuals", 
                               "Non-traditional Parenthood"))
seqdplot(match_seq, group = om5_match, xt = 15:60)

# 4 clusters:
om4_match <- cutree(match_clust, k = 4)
om4_match <- factor(om4_match, levels = c(4, 3, 2, 1),
                    labels = c("Early Family Formation", "Later Family Formation", 
                               "Childless Marriage", "Non-traditional Family Form"))
seqdplot(match_seq, group = om4_match, xt = 15:60)


summary(aov(as.numeric(om4_match) ~ as.factor(as.numeric(matches$group)))) 
# for pairwise 1:1, not significant for 4 and 5
# for propensity 1:1, not significant for 4 and 5
# for mahalanobis 1:1, not significant for 5 but for 4


# Regression analysis

# combine cluster membership

matches <- 
  matches %>%
  data.frame(., clust4 = om4_match, clust5 = om5_match) %>%
  mutate(gender = relevel(gender, ref = "Male"),
         country = relevel(country, ref = "France"),
         group = as.factor(as.numeric(group)),
         group = relevel(group, ref = "0"))

# rescaling? use resclae() by arm

# on 4 clusters
mult4 <- multinom(clust4 ~ gender + country + WW_age + group, 
                  data = matches)
stargazer(mult4, mult41, no.space = T, 
          covariate.labels = c("Female", "Austria", "Belgium", "Czech Republic", "Germany",
                               "Italy", "The Netherlands", "Poland", "Switzerland", "Age at WWII", "Treatment", "Constant"))

mult41 <- multinom(clust4 ~ gender + country + WW_age + evac + pow + camp, 
                  data = matches)
stargazer(mult4, no.space = T, 
          covariate.labels = c("Female", "Austria", "Belgium", "Czech Republic", "Germany",
                               "Italy", "The Netherlands", "Poland", "Switzerland", "Age at WWII", 
                               "Evacuation", "POW", "Camp", "Constant"))
# for pair 1:1, no group effect but some small exp effects
# propenstiy 1:1, small group effect plus exp effects
# for mahalanobis 1:1, very small group but some small exp effects

mult5 <- multinom(clust5 ~ gender + country + WW_age + group, 
                  data = matches)
stargazer(mult5, mult51, no.space=TRUE,
          covariate.labels = c("Female", "Austria", "Belgium", "Czech Republic", "Germany",
                               "Italy", "The Netherlands", "Poland", "Switzerland", "Age at WWII", "Treatment", "Constant"))

mult51 <- multinom(clust5 ~ country + gender + WW_age + evac + pow + camp + Wage * evac,
                       data = matches)
stargazer(mult5, no.space=TRUE,
          covariate.labels = c("Female", "Austria", "Belgium", "Czech Republic", "Germany",
                               "Italy", "The Netherlands", "Poland", "Switzerland", "Age at WWII", 
                               "Evacuation", "POW", "Camp", "Constant"))
# for pair 1:1, no group effect, some small evac effect
# for propensity 1:1, small group effect, some small evac and camp effect
# for mahalanobis 1:1, no group effect, evac effect

# z <- summary(res)$coefficients / summary(res)$standard.errors
# p <- (1 - pnorm(abs(z), 0, 1)) * 2
#stargazer(mult4, no.space=TRUE, 



#########
# i want to find all pairs where they are in the same/not in the same cluster

matches <- ddply(matches, .(propensity), transform, match = clust5[1] == clust5[2])

# this somewhat complicated method ensure differentiating in which cluster the treatment
# and which cluster the control group is

# filter to find the set of differing pairs
diff_pairs <-
  matches %>%
  subset(., match == F) %>%
  select(., group, propensity, clust5) 

# merging the data based on group ensures the differentiation between treatment and control
diff_pairs <- merge(diff_pairs[diff_pairs$group == 1,], diff_pairs[diff_pairs$group == 0,],
                    by = "propensity")
table(diff_pairs$clust5.x, diff_pairs$clust5.y)

# 526 pairs are in the same cluster
same_pairs <- 
  matches %>%
  subset(., match == T) %>%
  group_by(propensity) %>%
  arrange(clust5) %>%
  # show in which groups they are the same
  summarize(pair = paste(clust5, collapse = ""))
table(same_pairs$pair)


p_match <- glm(as.numeric(match) ~ gender + country + marage + clust5 + group, 
               data = matches, family = binomial(link = "logit"))


# plotting characteristic differences
library(ggplot2)

# country differences
ggplot(matches, aes(x = country)) + 
  geom_bar(aes(group = match, fill = match, y = ..prop..), position = "dodge") +
  scale_fill_manual(values = custom_col[c(2, 4)], labels = c("in different cluster", "in same cluster")) +
  labs(x = "Country", y = "Proportion", fill = "Treatment-control pair") + 
  theme_bw()
  
# gender differences
ggplot(matches, aes(x = gender, y = ..prop.., fill = match, group = match)) +
  geom_bar(position = "dodge") + 
  scale_fill_manual(values = custom_col[c(2, 4)], labels = c("in different cluster", "in same cluster")) +
  labs(x = "Gender", y = "Proportion", fill = "Treatment-control pair") +
  theme_bw()


# age diff
ggplot(matches, aes(x = WW_age, y = ..prop.., fill = match, group = match)) +
  geom_bar(position="identity", alpha = 0.5, width = 1) +
  scale_fill_manual(values = custom_col[c(2, 4)], labels = c("in different cluster", "in same cluster")) +
  labs(x = "Age at WWII", y = "Proportion", fill = "Treatment-control pair") +
  theme_bw()
  
matches$Wage <- ifelse(matches$WW_age < 18, 0, 1)
