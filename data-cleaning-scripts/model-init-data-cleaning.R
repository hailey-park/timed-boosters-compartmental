###################################################################################################
#Title: Cleaning Time-Since Data and Vaccine Data Used for Model Initialization
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
#Hospitalization data from this website (https://data.cdc.gov/Public-Health-Surveillance/Weekly-Rates-of-Laboratory-Confirmed-COVID-19-Hosp/6jg4-xsqq/about_data)
#Vaccine data from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Age-and-Sex-Trends-in-the-Uni/5i5k-6cmh/about_data)

hospitalizations <- read.csv("data/Weekly_Rates_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20240701.csv")
vaccines_cdc_age <- read.csv("data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional_20240624.csv")
combined_time_since_by_age <- read.csv("data/clean-data/combined_time_since_distributions_by_age.csv")[,-1]


#Total pop by age
total_pop_by_age <- vaccines_cdc_age %>% 
  filter(Location %in% c("CA"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years")) %>% group_by(Date, age_group) %>% summarise(total_pop = sum(census, na.rm=TRUE)) %>%
  group_by(age_group) %>% summarise(total_pop = mean(total_pop)) %>% mutate(pop_prop = total_pop/sum(total_pop))

write.csv(total_pop_by_age, "data/clean-data/total_pop_by_age.csv")[,-1]

#Severe cases (hospitalizations) at model initialization (Jan 1, 2024)
severe_cases_init <- hospitalizations %>% filter(State == "California", X_WeekendDate %in% c("2022-01-01", "2021-12-25", "2021-12-18"), AgeCategory_Legend %in% c("0-17 years (Children)", "18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                                 Race_Label == "All", Sex_Label == "All") %>% mutate(age_group = ifelse(AgeCategory_Legend == "0-17 years (Children)", "0-17 years", AgeCategory_Legend)) %>%
  select(X_WeekendDate, age_group, WeeklyRate) %>% group_by(age_group) %>% summarise(avg_inc = mean(WeeklyRate))#%>% dcast(age_group ~ X_WeekendDate) 

severe_cases_total_init <- merge(severe_cases_init, total_pop_by_age, by = c("age_group")) %>%
  mutate(severe_inf_total = avg_inc/100000 * total_pop) %>% select(age_group, severe_inf_total)

write.csv(severe_cases_total_init, "data/clean-data/severe_cases_total_init.csv")[,-1]

