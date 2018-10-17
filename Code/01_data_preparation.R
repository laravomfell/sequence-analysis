# This  file creates the original WWII survivor sequence dataset and writes it as a .csv
# It has two parts: 1. identifying the IDs of the survivors
#                   2. creating the family formation sequences from the Wave3 job episodes panel


# --------------
# Preliminaries
# --------------

# define path of the files
path <- "C:/Users/Lara/Dropbox/Uni/Sequence Analysis/Data/"

# load necessary libraries
library(dplyr)
library(foreign)

# -----------------------------
#
# 1. Identifying the survivors
#
# -----------------------------


# Identify the IDs of people which were directly affected by WWII
# by using AC002, Question 2 on the Accommodation Section

# load accommodation history file
accommodation <- read.dta(paste0(path, "Wave3_data/sharew3_rel6-0-0_ac.dta"))

# find people related to war. Here, non-eligible people (eg born after 1945 are not yet excluded)
war_subset <- accommodation %>%
            # evacuated during the war
  subset(sl_ac002d3 == "Selected" |
            # Prisoner of War
            sl_ac002d4 == "Selected" | 
            # concentration camp
            sl_ac002d7 == "Selected") %>%
  select(mergeid, sl_ac002d3, sl_ac002d4, sl_ac002d7) %>%
  # recoding the information
  mutate(evac = ifelse(sl_ac002d3 == "Selected", 1, 0),
         pow  = ifelse(sl_ac002d4 == "Selected", 1, 0),
         camp = ifelse(sl_ac002d7 == "Selected", 1, 0)) %>%
  # keep only relevant columns
  select(mergeid, evac, pow, camp)

# removing unneeded file
rm(accommodation); gc()

# ------------------------
#
# 2. create the sequences
#
# ------------------------

# read job episodes file containing the family formation data
jobs <- read.dta(paste0(path, "Wave3_data/sharewX_rel6-0-0_gv_job_episodes_panel.dta"))

# define the nine relevant countries
countries <- c("Austria", "Germany", "Netherlands", "Italy", "France", 
               "Switzerland", "Belgium", "Czech Republic", "Poland")

# create the sequence dataset
formation <- jobs %>% 
  # keep only entries of subgroup
  merge(war_subset) %>%
  # keep only participants born before 1944 
  # keep only the sequence information for ages 15-60
  # keep only participants born in the nine countries
  subset(yrbirth < 1945 & 
         age > 14 & 
         age < 61 & 
         country %in% countries) %>%
  # code family states: 10 = single, 20 = cohabitation, 30 = married
  mutate(state = ifelse(withpartner == 0 & married == 0, 10,
                        ifelse(withpartner == 1 & married == 0, 20,
                               ifelse(married == 1, 30, NA)))) %>%
  # add +1 if with children
  mutate(state = ifelse(nchildren > 0, state + 1, state)) %>%
  # drop variables not needed
  select(mergeid, hhid3, yrbirth, gender, age, year, country, state, evac, pow, camp) %>%
  # sort by ID and age to ensure correct ordering of sequences
  arrange(mergeid, age) %>%
  # recast data from long to wide format
  reshape(idvar = c("mergeid", "hhid3", "yrbirth", "gender", "country", "evac", "pow", "camp"), 
                 drop = "year", direction = "wide", timevar = "age") %>%
  # create WW_age variable and pre-allocate marriage age variable
  mutate(WW_age = 1944 - yrbirth,
         marage = NA)

# loop to find marriage ages
for (i in 15:60){
  # ugly fix to ensure correct looping through columns/ages
  col <- i - (15 - grep("state.", colnames(formation))[1])
  formation$marage <- ifelse(is.na(formation$marage) & (formation[,col] == 30 | formation[,col] == 31), i, formation$marage)
}

# clean up
rm(jobs, WW_subset, countries); gc()

# write file as .csv
write.csv(formation, file = paste0(path, "formation.csv"), row.names = F)