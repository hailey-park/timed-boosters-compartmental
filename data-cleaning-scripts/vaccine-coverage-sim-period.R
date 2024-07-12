###################################################################################################
#Title: Cleaning Vaccine Coverage Data Used
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
#Vaccine data through May 10, 2023 from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Age-and-Sex-Trends-in-the-Uni/5i5k-6cmh/about_data)
#Vaccine data from June 14 - October 11 2023 from this website (https://data.cdc.gov/Vaccinations/COVID-19-Vaccines-Up-to-Date-Status/9b5z-wnve/about_data)
#Vaccine data from July 1 2023 - March 31 2024 from this website (https://data.cdc.gov/Vaccinations/Monthly-Cumulative-Number-and-Percent-of-Persons-W/vugp-mqip/about_data)

vaccines_may2023 <- read.csv("data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional_20240624.csv")
vaccines_oct2023 <- read.csv("data/COVID-19_Vaccines_Up_to_Date_Status_20240702.csv")
vaccine_mar2024 <- read.csv("data/Monthly_Cumulative_Number_and_Percent_of_Persons_Who_Received_1__updated_2023-24_COVID-19_Vaccination_Doses_by_Age_Group_and_Jurisdiction__United_States_20240702.csv")


#Vaccine doses by week through May 2023
doses_by_week_may2023 <- vaccines_may2023 %>% 
  filter(Location %in% c("CA"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y%H"),
         week = floor_date(date, unit = "week")) %>%
  filter(date >= as.Date("2021-12-25") & date <= as.Date('2024-06-01')) %>%
  group_by(age_group, date, week) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                                total_primary_series = sum(Series_Complete_Yes, na.rm=TRUE),
                                                total_boosted_1st = sum(Booster_Doses, na.rm=TRUE),
                                                total_boosted_2nd = sum(Second_Booster, na.rm=TRUE)) %>%
  group_by(age_group, week) %>% summarise(total_pop = mean(total_pop, na.rm=TRUE),
                                          total_primary_series = max(total_primary_series),
                                          total_boosted_1st = max(total_boosted_1st),
                                          total_boosted_2nd = max(total_boosted_2nd)) %>%
  group_by(age_group) %>% mutate(primary_series_today = c(total_primary_series[1], (total_primary_series - lag(total_primary_series, 1))[-1]),
                                 boosted_1st_today = c(total_boosted_1st[1], (total_boosted_1st - lag(total_boosted_1st, 1))[-1]),
                                 boosted_2nd_today = c(total_boosted_2nd[1], (total_boosted_2nd - lag(total_boosted_2nd, 1))[-1]),
                                 primary_series_rate = primary_series_today/total_pop,
                                 boosted_rate = (boosted_1st_today + boosted_2nd_today)/total_primary_series,
                                 week_num = as.numeric(difftime(week, as.Date("2021-12-26"), units = "weeks"))) %>%
  filter(week != as.Date("2021-12-19"))

write.csv(doses_by_week_may2023 %>% select(age_group, primary_series_rate, week, week_num), "data/clean-data/primary_vax_rate_over_time.csv")
write.csv(doses_by_week_may2023 %>% select(age_group, boosted_rate, week, week_num), "data/clean-data/booster_vax_rate_over_time.csv")


#Vaccine doses by month from June 2023 through October 2023
doses_by_month_oct2023 <- vaccines_oct2023 %>% 
  filter(Location %in% c("CA"), 
         Demographic_Category %in% c("Ages_75+_yrs", "Ages_65-74_yrs", "Ages_50-64_yrs", "Ages_25-49_yrs", "Ages_18-24_yrs", "Ages_12-17_yrs", "Ages_5-11_yrs", "Ages_<5yrs")) %>% 
  mutate(age_group = case_when(Demographic_Category %in% c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs") ~ "0-17 years",
                               Demographic_Category %in% c("Ages_18-24_yrs", "Ages_25-49_yrs") ~ "18-49 years",
                               Demographic_Category == "Ages_50-64_yrs" ~ "50-64 years",
                               Demographic_Category == "Ages_65-74_yrs" ~ "65-74 years",
                               Demographic_Category == "Ages_75+_yrs" ~ "75+ years"),
         date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(date >= as.Date('2022-01-01') & date <= as.Date('2024-06-01')) %>%
  group_by(age_group, date) %>% summarise(total_pop = sum(census, na.rm=TRUE),
                                                total_up_to_date = sum(Up_to_date, na.rm=TRUE)) %>%
  group_by(age_group) %>% mutate(up_to_date_today = c((total_up_to_date - lag(total_up_to_date, 1))[2], (total_up_to_date - lag(total_up_to_date, 1))[-1]))

#Vaccine doses by month from Jul 2023 through March 2024
doses_by_month_mar2024 <- vaccine_mar2024 %>% filter(Jurisdiction == "California", Age_group_label %in% c("6 months to 17 years", "18-49 years", "50-64 years", "65+ years")) %>%
  select(Month, Numerator, Population, Age_group_label)


