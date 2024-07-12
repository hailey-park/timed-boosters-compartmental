###################################################################################################
#Title: Cleaning Case Data and Vaccine Data Used for 'Time Since Last' Estimation
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


#Read in data
#Case data from this website (https://data.chhs.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state)
#Vaccine data from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Age-and-Sex-Trends-in-the-Uni/5i5k-6cmh/about_data)

cases <- read.csv("data/covid19casesdemographics.csv")
vaccines_cdc_age <- read.csv("data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional_20240624.csv")

#Cases by Week
cases_by_week <- cases %>% mutate(report_date = as.Date(report_date)) %>% 
  filter(demographic_category == "Age Group", report_date >= as.Date('2020-03-01') & report_date <= as.Date('2021-12-31'), demographic_value %in% c("0-17", "18-49", "50-64", "65+")) %>%
  group_by(demographic_value) %>% mutate(lag_cases = c(0, diff(total_cases, 1)), lag_cases = if_else(lag_cases < 0, 0, lag_cases)) %>% 
  mutate(week = floor_date(report_date, unit = "week")) %>% group_by(week, demographic_value) %>% summarise(total_cases = sum(lag_cases)) %>% group_by(demographic_value) %>%
  mutate(perc_cases = total_cases/sum(total_cases),
         demographic_value = case_when(demographic_value == "65+" ~ "65-74 years",
                                       demographic_value == "50-64" ~ "50-64 years",
                                       demographic_value == "18-49" ~ "18-49 years",
                                       demographic_value == "0-17" ~ "0-17 years"))

cases_75plus <- cases_by_week %>% filter(demographic_value == "65-74 years") %>% mutate(demographic_value = "75+ years")

cases_by_week <- as.data.frame(rbind(cases_by_week, cases_75plus))


#Vaccine Doses by Week
#NOTE: We will be shortening the distribution for fully vaccinated doses to include individuals
#      who only recently completed their primary-series (excluding primary series before March
#      22, 2021 -- 6 months before from when first boosters doses are administered)

doses_by_week <- vaccines_cdc_age %>% 
  filter(Location %in% c("CA"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y%H")) %>%
  filter(date >= as.Date('2020-12-15') & date <= as.Date('2021-12-31')) %>%
  group_by(age_group, date) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                          total_primary_series = sum(Series_Complete_Yes, na.rm=TRUE),
                                          total_boosted = sum(Booster_Doses, na.rm=TRUE)) %>%
  group_by(age_group) %>% mutate(primary_series_today = c(total_primary_series[1], (total_primary_series - lag(total_primary_series, 1))[-1]),
                                 boosted_today = c(total_boosted[1], (total_boosted - lag(total_boosted, 1))[-1]))

fully_vax_doses_by_week <- doses_by_week %>% #filter(date >= as.Date("2021-03-22")) %>%
  mutate(week = floor_date(as.Date(date), unit = "week")) %>% group_by(week, age_group) %>% summarise(total_vaccinated = sum(primary_series_today)) %>%
     group_by(age_group) %>% mutate(perc_doses = total_vaccinated/sum(total_vaccinated))

booster_doses_by_week <- doses_by_week %>% filter(date >= as.Date("2021-09-22")) %>%
  mutate(week = floor_date(as.Date(date), unit = "week")) %>% group_by(week, age_group) %>% summarise(total_vaccinated = sum(boosted_today)) %>%
  group_by(age_group) %>% mutate(perc_doses = total_vaccinated/sum(total_vaccinated))

#write as .csv
write.csv(cases_by_week, "data/clean-data/cases_by_week.csv")
write.csv(booster_doses_by_week, "data/clean-data/booster_doses_by_week.csv")
write.csv(fully_vax_doses_by_week, "data/clean-data/fully_vax_doses_by_week.csv")

#Plot coverages over time
ggplot(data = doses_by_week, aes(x = date, y = total_primary_series/total_pop, color = age_group)) +
  geom_line() +
  xlab("Time") +
  ylab("Proportion") +
  ggtitle("Proportion of Completed Primary Series Over Time")
  

#Getting age-specific vaccine coverage
vaccine_coverage_cdc_age <- vaccines_cdc_age %>% filter(Date == "12/31/2021 12:00:00 AM", Location %in% c("US", "CA"), Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs") ) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years")) %>%
  group_by(age_group, Location) %>%
  summarise(census = sum(census),
            Administered_Dose1 = sum(Administered_Dose1, na.rm=TRUE),
            Series_Complete_Yes = sum(Series_Complete_Yes, na.rm=TRUE),
            Booster_Doses = sum(Booster_Doses, na.rm=TRUE)) %>% mutate(booster_coverage_eligible = Booster_Doses/Series_Complete_Yes,
                                                           booster_coverage = Booster_Doses/census,
                                                           fully_vacc_coverage = Series_Complete_Yes/census)


# vaccines <- read.csv("data/covid19vaccinesadministeredbydemographics-booster-archived.csv")
# vaccines_cdc <- read.csv("data/COVID-19_Vaccinations_in_the_United_States_Jurisdiction_20240624.csv")
  # full_vax_doses_by_week <- vaccines %>% mutate(week = floor_date(as.Date(administered_date), unit = "week"), administered_date = as.Date(administered_date)) %>% 
  #   filter((administered_date >= as.Date('2020-12-15') & administered_date <= as.Date('2021-12-01')), demographic_category == "Age Group", demographic_value %in% c("Under 5", "5-11", "12-17", "18-49", "50-64", "65+")) %>%
  #   mutate(age_group = if_else(demographic_value %in% c("Under 5", "5-11", "12-17"), "0-17", demographic_value)) %>% group_by(week, age_group) %>% summarise(total_vaccinated = sum(fully_vaccinated)) %>%
  #   group_by(age_group) %>% mutate(perc_doses = total_vaccinated/sum(total_vaccinated))
  #
  # booster_doses_by_week <- vaccines %>% mutate(week = floor_date(as.Date(administered_date), unit = "week"), administered_date = as.Date(administered_date)) %>% 
  #   filter((administered_date >= as.Date('2021-09-22') & administered_date <= as.Date('2021-12-01')), demographic_category == "Age Group", demographic_value %in% c("Under 5", "5-11", "12-17", "18-49", "50-64", "65+")) %>%
  #   mutate(age_group = if_else(demographic_value %in% c("Under 5", "5-11", "12-17"), "0-17", demographic_value)) %>% group_by(week, age_group) %>% summarise(total_boosted = sum(booster_recip_count)) %>%
  #   group_by(age_group) %>% mutate(perc_doses = total_boosted/sum(total_boosted))
  
# booster_vaccine_coverage_cdc <- vaccines_cdc %>% filter(Date == "12/31/2021", Location %in% c("US", "CA")) %>% 
#   dplyr::select(Date, Location, Additional_Doses_Vax_Pct, Additional_Doses_18Plus_Vax_Pct, Additional_Doses_50Plus_Vax_Pct, Additional_Doses_65Plus_Vax_Pct)

# full_vax_vaccine_coverage <- vaccines %>% mutate(week = floor_date(as.Date(administered_date), unit = "week"), administered_date = as.Date(administered_date)) %>% 
#   filter((administered_date == as.Date('2021-12-01')), demographic_category == "Age Group", demographic_value %in% c("Under 5", "5-11", "12-17", "18-49", "50-64", "65+")) %>%
#   mutate(age_group = if_else(demographic_value %in% c("Under 5", "5-11", "12-17"), "0-17", demographic_value)) %>% group_by(age_group) %>% summarise(cumulative_fully_vaccinated = sum(cumulative_fully_vaccinated)) %>%
#   mutate(total_pop = 39240000 * c(.2176, .4419, .1825, .158), coverage = cumulative_fully_vaccinated/total_pop)
# 
# booster_vaccine_coverage <- vaccines %>% mutate(week = floor_date(as.Date(administered_date), unit = "week"), administered_date = as.Date(administered_date)) %>% 
#   filter((administered_date == as.Date('2022-01-31')), demographic_category == "Age Group", demographic_value %in% c("Under 5", "5-11", "12-17", "18-49", "50-64", "65+")) %>%
#   mutate(age_group = if_else(demographic_value %in% c("Under 5", "5-11", "12-17"), "0-17", demographic_value)) %>% group_by(age_group) %>% summarise(cumulative_boosted = sum(cumulative_booster_recip_count),
#                                                                                                                                                      boosted_eligible = sum(booster_eligible_population),
#                                                                                                                                                      cumulative_fully_vaccinated = sum(cumulative_fully_vaccinated)) %>%
#   mutate(total_pop = 39240000 * c(.2176, .4419, .1825, .158), coverage = cumulative_boosted/total_pop,
#          coverage_eligible = cumulative_boosted/boosted_eligible,
#          coverage_eligible_alt = cumulative_boosted/cumulative_fully_vaccinated)
# 
