###################################################################################################
#Title: Cleaning Incidence Data for Model Calibration
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

hospitalizations <- read.csv("data/Weekly_Rates_of_Laboratory-Confirmed_COVID-19_Hospitalizations_from_the_COVID-NET_Surveillance_System_20240701.csv")

#Severe cases (hospitalizations) at model initialization (Jan 1, 2024)
severe_cases_inc <- hospitalizations %>% mutate(age_group = ifelse(AgeCategory_Legend == "0-17 years (Children)", "0-17 years", AgeCategory_Legend),
                                                     week = as.Date(X_WeekendDate)) %>% 
  filter(State == "California", AgeCategory_Legend %in% c("0-17 years (Children)", "18-49 years", "50-64 years", "65-74 years", "75+ years"), Race_Label == "All", Sex_Label == "All",
         week >= as.Date("2022-01-01") & week <= as.Date("2023-05-10")) %>%
  select(week, age_group, WeeklyRate) %>% dcast(week ~ age_group) 

write.csv(severe_cases_inc, "data/clean-data/severe_cases_inc.csv")[,-1]
