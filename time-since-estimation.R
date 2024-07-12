###################################################################################################
#Title: 'Time Since Last' Estimation
#Author: Hailey Park
#Date: June 19, 2024
###################################################################################################

rm(list=ls())

setwd(here::here())

#Load libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)

#Read data
cases_by_week <- read.csv("data/clean-data/cases_by_week.csv")[,-1]
booster_doses_by_week <- read.csv("data/clean-data/booster_doses_by_week.csv")[,-1]
fully_vax_doses_by_week <- read.csv("data/clean-data/fully_vax_doses_by_week.csv")[,-1]

#'Time Since Last' Estimation -- Initializing model at December 31, 2021
#NOTE: Estimation run separately for each age group because of differences in vaccine coverage and seroprevalence

#Clean data
booster_doses_by_week$week <- as.character(booster_doses_by_week$week)
fully_vax_doses_by_week$week <- as.character(fully_vax_doses_by_week$week)
cases_by_week$week <- as.character(cases_by_week$week)

#Calculate time since last vaccine dose or infection
add.weeks= function(date,n) {seq(date, by = paste (n, "weeks"), length = 2)[2]}

time_since_last <- function(df) {
  
  #Get age-specific dose/case distributions
  booster_doses_by_week <- booster_doses_by_week %>% filter(age_group == df$age_group[1])
  fully_vax_doses_by_week <- fully_vax_doses_by_week %>% filter(age_group == df$age_group[1])
  cases_by_week <- cases_by_week %>% filter(demographic_value == df$age_group[1])
  
  #Calculate time since last dose and time since last infection
  set.seed(88)
  last_dose_and_inf <- df %>% mutate(time_since_last_dose = ifelse(prior_vacc == 'fullvax', sample(as.character(fully_vax_doses_by_week$week), size = sum(prior_vacc == 'fullvax'), prob = fully_vax_doses_by_week$perc_doses, replace = TRUE),
                                                                   ifelse(prior_vacc == 'booster', sample(as.character(booster_doses_by_week$week), size = sum(prior_vacc == 'booster'), prob = booster_doses_by_week$perc_doses, replace = TRUE),
                                                                          NA)),
                                     time_since_last_inf = ifelse(prior_inf == 1,
                                                                  sample(as.character(cases_by_week$week), 
                                                                         size = sum(prior_inf == 1),
                                                                         prob = cases_by_week$perc_cases,
                                                                         replace = TRUE),
                                                                  NA),
                                     time_since_last_dose_inf = pmax(as.Date(time_since_last_dose), as.Date(time_since_last_inf), na.rm =  TRUE),
                                     immunity_type = case_when(prior_inf == 1 ~ "hybrid/infection-only",
                                                               prior_vacc != 'unvax' & prior_inf == 0 ~ "vaccine-only",
                                                               prior_vacc == 'unvax' & prior_inf == 0 ~ "immune-naive",
                                                               TRUE ~ NA),
                                     time_since = difftime(as.Date('2021-12-31'), time_since_last_dose_inf, units = "days"),
                                     time_since_inf = difftime(as.Date('2021-12-31'), time_since_last_inf, units = "days"),
                                     waning_group = case_when(time_since < 90 ~ "0-3 months",
                                                              time_since >= 90 & time_since < 182 ~ "3-6 months",
                                                              time_since >= 182 & time_since < 365 ~ "6-12 months",
                                                              time_since >= 365 & time_since < 730 ~ "12+ months",
                                                              TRUE ~ NA),
                                     immunity_type = ifelse(immunity_type == "hybrid/infection-only" & time_since_inf < 90, "recovered", immunity_type),
                                     waning_group = ifelse(immunity_type == "hybrid/infection-only" & time_since_inf < 90, NA, waning_group)) 

  return(last_dose_and_inf)
}

#Create matrices for each age group
#NOTE: Seroprevalence estimates are age-specific and at the state-level (CA), using Dec 2021.
#.     https://covid.cdc.gov/covid-data-tracker/#national-lab
#.     Prior vaccine estimates are age-specific and at the state-level (CA), using cumulative fully vaccinated estimates at December 31, 2021
#.     https://data.chhs.ca.gov/dataset/vaccine-progress-dashboard

set.seed(88)
age_0_17_cal <- data.frame(individual = c(1:1000000),
                           age_group = '0-17 years',
                           prior_vacc = as.character(sample(c('unvax', 'fullvax', 'booster'), 1000000, prob = c((1 - 0.309), (0.309 - 0.001), 0.001), replace = TRUE)),
                           prior_inf = rbinom(1000000, 1, 0.417))
set.seed(88)
age_18_49_cal <- data.frame(individual = c(1:1000000),
                            age_group = '18-49 years',
                            prior_vacc = as.character(sample(c('unvax', 'fullvax', 'booster'), 1000000, prob = c((1 - 0.735), (0.735 - 0.224), 0.224), replace = TRUE)),
                            prior_inf = rbinom(1000000, 1, 0.344))
set.seed(88)
age_50_64_cal <- data.frame(individual = c(1:1000000),
                            age_group = '50-64 years',
                            prior_vacc = as.character(sample(c('unvax', 'fullvax', 'booster'), 1000000, prob = c((1 - 0.843), (0.843 - 0.393), 0.393), replace = TRUE)),
                            prior_inf = rbinom(1000000, 1, 0.318))

set.seed(88)
age_65_74_cal <- data.frame(individual = c(1:1000000),
                            age_group = '65-74 years',
                            prior_vacc = as.character(sample(c('unvax', 'fullvax', 'booster'), 1000000, prob = c((1 - 0.912), (0.912 - 0.579), 0.579), replace = TRUE)),
                            prior_inf = rbinom(1000000, 1, 0.196))

set.seed(88)
age_75plus_cal <- data.frame(individual = c(1:1000000),
                             age_group = '75+ years',
                             prior_vacc = as.character(sample(c('unvax', 'fullvax', 'booster'), 1000000, prob = c((1 - 0.841), (0.841 - 0.583), 0.583), replace = TRUE)),
                             prior_inf = rbinom(1000000, 1, 0.196))


set.seed(88)
time_since_results <- list(age_0_17_cal, age_18_49_cal, age_50_64_cal, age_65_74_cal, age_75plus_cal) %>%
  lapply(time_since_last)

inspection <- time_since_results[[5]] %>% group_by(age_group, immunity_type, waning_group) %>% summarise(proportion = n()/1000000)



#population estimates from CDC vaccine data of CA
combine_results <- time_since_results %>%
  lapply(function(df){df %>% group_by(age_group, immunity_type, waning_group) %>% summarise(proportion = n()/1000000)})
combined <- bind_rows(combine_results) %>% mutate(num_people = round(case_when(age_group == "0-17 years" ~ 8894641 * proportion,
                                                                         age_group == "18-49 years" ~ 17528506 * proportion,
                                                                         age_group == "50-64 years" ~ 7250961 * proportion,
                                                                         age_group == "65-74 years" ~ 3386670 * proportion,
                                                                         age_group == "75+ years" ~ 2451445 * proportion)))

write.csv(combined, "data/clean-data/combined_time_since_distributions_by_age.csv")
############################################################
#PLOTS
inspection <- time_since_results[[1]] %>% select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_dose_inf) %>%
  group_by(prior_vacc, time_since_last_dose) %>% summarise(total = n()) %>% filter(prior_vacc != 'unvax')

observed_data <- as.data.frame(rbind(fully_vax_doses_by_week %>% mutate(prior_vacc = "observed_fullvax"),booster_doses_by_week %>% mutate(prior_vacc = "observed_boosted"))) %>% 
  select(prior_vacc, week, age_group, perc_doses) %>% filter(age_group == "0-17 years")

#plot simulated prior dose
ggplot(data = inspection, aes(x = as.Date(time_since_last_dose), y = total, color = prior_vacc)) +
  geom_line() +
  ylab("Total COVID-19 Doses") +
  xlab("Weeks") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2020-12-15"), as.Date("2021-12-31")),
               breaks = "1 months") +
  ylim(0, 120000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Simulated Time-Since Last Vaccination\nAge Group: 0-17 years")

#plot simulated prior dose with observed data
ggplot() +
  geom_line(data = inspection, aes(x = as.Date(time_since_last_dose), y = total, color = prior_vacc)) +
  geom_line(data = observed_data, aes(x = as.Date(week), y = perc_doses * 500, color = prior_vacc)) +
  xlab("Weeks") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2020-12-15"), as.Date("2021-12-31")),
               breaks = "1 months") +
  ylim(0, 300000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Simulated Time-Since Last Vaccination\nAge Group: 0-17 years") +
  scale_y_continuous(
    # Features of the first axis
    name = "Total COVID-19 Doses",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~./500, name="Observed Proportion\n(Denominator by vaccine status)")
  )



inspection <- time_since_results[[5]] %>% select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_dose_inf) %>%
  group_by(time_since_last_inf) %>% summarise(total = n()) %>% filter(!is.na(time_since_last_inf))

observed_data <- cases_by_week %>% filter(demographic_value == "75+ years")
  
## plot simulated prior inf
# ggplot(data = inspection, aes(x = as.Date(time_since_last_inf), y = total)) +
#   ylab("Total COVID-19 Infections") +
#   xlab("Time") +
#   scale_x_date(date_labels = ("%b-%Y"),
#                limits = c(as.Date("2020-03-01"), as.Date("2021-12-31")),
#                breaks = "1 months") +
#   ylim(0, 10000)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("Simulated Time-Since Last Infection\nAge Group: 75+ years")

#plot simulated prior inf with observed data
ggplot() +
  geom_line(data = inspection, aes(x = as.Date(time_since_last_inf), y = total, color = "Simulated")) +
  geom_line(data = observed_data, aes(x = as.Date(week), y = perc_cases * 200000, color = "Observed")) +
  
  xlab("Time") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2020-03-01"), as.Date("2021-12-31")),
               breaks = "1 months") +
  ylim(0, 10000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Simulated Time-Since Last Infection\nAge Group: 75+ years") +
  scale_y_continuous(
    # Features of the first axis
    name = "Total Simulated COVID-19 Infections",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans=~./200000, name="Observed Proportion")
  )


inspection <- time_since_results[[5]] %>% select(prior_vacc, prior_inf, time_since_last_dose, time_since_last_inf, time_since_last_dose_inf) %>%
  group_by(time_since_last_dose_inf) %>% summarise(total = n()) %>% filter(!is.na(time_since_last_dose_inf))

#plot simulated time since
ggplot(data = inspection, aes(x = as.Date(time_since_last_dose_inf), y = total)) +
  geom_line() +
  ylab("Total COVID-19 Immune Events") +
  xlab("Time") +
  scale_x_date(date_labels = ("%b-%Y"),
               limits = c(as.Date("2020-03-01"), as.Date("2021-12-31")),
               breaks = "1 months") +
  ylim(0, 110000)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Simulated Time-Since Last Immune Event\nAge Group: 75+ years")


