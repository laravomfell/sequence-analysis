# This file does the following:
# 1. extract the IDs of potential control cases from the full dataset
# 2. consider a range of different treatment-control matchings and check 
#    them for data balance
# 3. pick the best matching
# 4. create the final sequence set and save it as a csv


# Preliminaries

# define path of the files
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# load necessary libraries
library(dplyr)     # data transformations
library(foreign)   # to load stata file
library(RItools)   # for balance checks
library(optmatch)  # for matching


# ---------------------------
#                           
# 1. FIND POTENTIAL MATCHES 
#                           
# ---------------------------


# First, load data on WWII survivors
formation <- read.csv(paste0(path, "formation.csv"))
# create "treatment" dummy for WW survivors
formation$group <- T

# read job episodes file containing the family formation data
jobs <- read.dta(paste0(path, "sharewX_rel6-0-0_gv_job_episodes_panel_stata/sharewX_rel6-0-0_gv_job_episodes_panel.dta"))

# define relevant countries
countries <- c("Austria", "Germany", "Netherlands", "Italy", "France", "Switzerland",
               "Belgium", "Czech Republic", "Poland")

# create dataframe with non-survivor information
non_WWII_ID <- jobs %>% 
  # delete non WWII participants
  filter(., !(mergeid %in% formation$mergeid)) %>%
  # remove all participants born after 1944 and in the non-relevant countries
  subset(., yrbirth < 1945 &
           country %in% countries &
           !duplicated(mergeid)) %>%
  # drop variables not needed
  select(mergeid, yrbirth, gender, country) %>%
  # create "treatment" dummy
  mutate(group = F,
         gender = droplevels(gender),
         country = droplevels(country))

# merge data, keep only columns relevant for matching
final <- rbind(non_WWII_ID, formation[, colnames(formation) %in% names(non_WWII_ID)])

# check for balancing of data
xBalance(group ~ . -mergeid, 
         data = final, 
         report = c("chisquare.test")) # not balanced


# -------------------------------
#
# 2. CREATE DIFFERENT MATCHINGS 
#
# -------------------------------

# start with simple pairwise matching but with two controls per case
pairs <- pairmatch(group ~ yrbirth + gender + country, 
                   controls = 1, data = final)

# propensity score matching
psm <- glm(group ~ yrbirth + gender + country, family = "binomial", data = final)
# check sufficient overlap
boxplot(psm) # some overlap

# create propensity score match
propensity <- pairmatch(psm, data = final, controls = 1)

# check difference between pairwise and propensity matching
all.equal(pairs, propensity, check.attributes = F) # not the same

# Use Mahalanobis distance
# create distance 
mahalanobis <- match_on(group ~ yrbirth + scores(psm), data = final)
# create matching
maha_pairs <- pairmatch(mahalanobis, caliper = 1, data = final, controls = 1)
# save up some memory by deleting distance matrix
rm(mahalanobis); gc()

# check differences
all.equal(maha_pairs, pairs, check.attributes = F) # not the same
all.equal(propensity, maha_pairs, check.attributes = F) # not the same

# Balance testing: plot differences between treatment and control
allbalance <- xBalance(group ~ . -mergeid, 
                       data = final, 
                       report = c("chisquare.test", "std.diffs"), 
                       strata = data.frame(Original = factor("none"), 
                                           Pairwise = pairs, 
                                           Propensity = propensity, 
                                           Mahalanobis = maha_pairs))

# --------------------------------------
# This part is only needed for plotting

library(ggplot2)
library(viridis)
library(tidyr)

# fortify as dataframe
balance_plot <- as.data.frame(RItools:::prepareXbalForPlot(allbalance))

# rename variables for plot
balance_plot$name <- factor(c("Birth Year", "Male", "Female", "Austria", 
                              "Germany", "Netherland", "Italy", "France",
                              "Switzerland", "Belgium", "Czech Republic", "Poland"))
# reorder dataframe
balance_plot <- gather(balance_plot, strata, values, -name)

# create plot
ggplot(balance_plot, aes(y = name, x = values, color = strata)) + 
  geom_vline(xintercept = 0) +
  geom_point(size = 3) +
  labs(x = "Standardized Differences", y = NULL, color = "Matchings") +
  scale_color_viridis(discrete = T, 
                      limits=c("Original", "Pairwise", "Propensity", "Mahalanobis")) +
  theme(legend.position = "bottom") +
  theme_bw() +
  scale_y_discrete(limits = rev(c("Birth Year", "Male", "Female", "Austria", "Germany", "Netherland", "Italy", "France",
                            "Switzerland", "Belgium", "Czech Republic", "Poland")))
# -----------------------------------------------------------------------------------


# -------------------
# 
# 3. Interpretation  
#
# -------------------

# pairwise is best fit, except for birth year
# next best: propensity score matching

#-----------------------------
#
# 4. create matching dataset
#
#------------------------------

match_ID <- 
  final %>%
  cbind(., propensity) %>%
  subset(., complete.cases(.) & group == F)

match_sequences <- 
  jobs %>%
  subset(., mergeid %in% match_ID$mergeid &
            yrbirth < 1945 & 
            age > 14 & 
            age < 61) %>%
  # code family states
  mutate(state = ifelse(withpartner == 0 & married == 0, 10,
                        ifelse(withpartner == 1 & married == 0, 20,
                               ifelse(married == 1, 30, NA)))) %>%
  # add +1 if with children
  mutate(state = ifelse(nchildren > 0, state + 1, state),
         # add non WWII surviving variables
         evac = 0,
         pow = 0,
         camp = 0, 
         group = F) %>%
  # drop variables not needed
  select(., mergeid, hhid3, yrbirth, gender, age, year, country, state, evac, pow, camp, group) %>%
  # sort by ID and age to ensure correct ordering of sequences
  arrange(., mergeid, age) %>%
  # recast data from long to wide format
  reshape(., idvar = c("mergeid", "hhid3", "yrbirth", "gender", "country", "evac", "pow", "camp", "group"), 
          drop = "year", direction = "wide", timevar = "age") %>%
  # create WW_age variable and marriage age variable
  mutate(WW_age = 1944 - yrbirth,
         marage = NA)

# loop to find marriage ages
for (i in 15:60){
  cols <- i - 5
  match_sequences$marage <- ifelse(is.na(match_sequences$marage) & 
                                   (match_sequences[,cols] == 30 | match_sequences[,cols] == 31), 
                                   i, match_sequences$marage)
}

matches <- rbind(match_sequences, formation)
matches$propensity <- as.character(propensity[!is.na(propensity)])

write.csv(matches, file = paste0(path, "matches.csv"), row.names = F)

# --------------------------------------------
#
# Plotting group differences between matches
# 
# --------------------------------------------

# define colors
custom_col <- viridis(n = 6, begin = 1, end = 0)

# plotting gender differences

# need helper table for correct proportions
helper <- as.data.frame(prop.table(table(matches$gender, matches$group), margin = 2))
names(helper) <- c("Gender", "Group", "Freq")

# plot gender differences
ggplot(helper, aes(x = Group, y = Freq, fill = Gender)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_col[c(2, 4)]) +
  scale_x_discrete(labels = c("Control", "Treatment")) +
  labs(y = "Proportion") +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
  theme_bw()

# country differences
ggplot(matches, aes(x = country)) + 
  geom_bar(aes(group = group, fill = group, y = ..prop..), position = "dodge") +
  scale_fill_manual(values = custom_col[c(2, 4)], labels = c("in different cluster", "in same cluster")) +
  labs(x = "Country", y = "Proportion", fill = "Treatment-control pair") + 
  theme_bw()

# age differences
ggplot(matches, aes(x = WW_age, y = ..prop.., fill = group, group = group)) +
  geom_bar(position="identity", alpha = 0.5, width = 1) +
  scale_fill_manual(values = custom_col[c(2, 4)], labels = c("Control", "Treatment")) +
  labs(x = "Age at WWII", y = "Proportion", fill = "Group") +
  theme_bw()