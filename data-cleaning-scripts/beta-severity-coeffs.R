###################################################################################################
#Title: Script for creating coefficients for beta and severity terms (match protection curves)
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

#Read in waning curve data
waning_curves_nonsevere <- read.csv("data/clean-data/combined_nonsevere_waning_predictions_weekly.csv")[,-1]
waning_curves_severe <- read.csv("data/clean-data/combined_severe_waning_predictions_weekly.csv")[,-1]

#Adjust waning curves with the model input parameters
adj_waning_curves_nonsevere <- waning_curves_nonsevere %>% filter(immunocompromised == 0, estimate == "mean") %>%
  mutate(waning_stage = case_when(weeks %in% c(1:12) ~ "0-3 months",
                                  weeks %in% c(13:25) ~ "3-6 months",
                                  weeks %in% c(26:51) ~ "6-12 months",
                                  weeks %in% c(52:104) ~ "12+ months")) %>% rowwise() %>%
  mutate(ve_pred = max(ve_pred, 0)) %>% group_by(age_group, prior_inf, waning_stage) %>% summarise(avg_ve_pred = mean(ve_pred)) %>% arrange(age_group)

adj_waning_curves_severe <- waning_curves_severe %>% filter(immunocompromised == 0, estimate == "mean") %>%
  mutate(waning_stage = case_when(weeks %in% c(1:12) ~ "0-3 months",
                                  weeks %in% c(13:25) ~ "3-6 months",
                                  weeks %in% c(26:51) ~ "6-12 months",
                                  weeks %in% c(52:104) ~ "12+ months")) %>% rowwise() %>%
  mutate(ve_pred = max(ve_pred, 0)) %>% group_by(age_group, prior_inf, waning_stage) %>% summarise(avg_ve_pred = mean(ve_pred)) %>% arrange(age_group)

#Calculate beta coefficients
beta_coeffs <- data.frame(age_group = c("0-17 years", "18-49 years", "50-64 years", "65-74 years", "75+ years"),
                    beta_0 = 1,              
                    beta_V_1 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 0, waning_stage == "0-3 months"))$avg_ve_pred),
                    beta_V_2 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 0, waning_stage == "3-6 months"))$avg_ve_pred), 
                    beta_V_3 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 0, waning_stage == "6-12 months"))$avg_ve_pred),
                    beta_V_4 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 0, waning_stage == "12+ months"))$avg_ve_pred),
                    beta_1 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 1, waning_stage == "0-3 months"))$avg_ve_pred),
                    beta_2 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 1, waning_stage == "3-6 months"))$avg_ve_pred),
                    beta_3 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 1, waning_stage == "6-12 months"))$avg_ve_pred),
                    beta_4 = 1 - ((adj_waning_curves_nonsevere %>% filter(prior_inf == 1, waning_stage == "12+ months"))$avg_ve_pred)) 

write.csv(beta_coeffs, "data/clean-data/beta_coeffs.csv")

beta_coefs_relabeled <- melt(beta_coeffs) %>% mutate(immunity_type = case_when(variable == "beta_0" ~ "immune-naive",
                                                                               grepl("V", variable) ~ "vaccine-only",
                                                                               TRUE ~ "hybrid/infection-only"),
                                                     waning_group = case_when(grepl("1", variable) ~ "0-3 months",
                                                                              grepl("2", variable) ~ "3-6 months",
                                                                              grepl("3", variable) ~ "6-12 months",
                                                                              grepl("4", variable) ~ "12+ months",
                                                                              TRUE ~ NA))
write.csv(beta_coefs_relabeled, "data/clean-data/beta_coeffs_relabeled.csv")


severity_coeffs <- data.frame(age_group = c("0-17 years", "18-49 years", "50-64 years", "65-74 years", "75+ years"),
                    severity_0 = 1,              
                    severity_V_1 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 0, waning_stage == "0-3 months"))$avg_ve_pred))/beta_coeffs$beta_V_1,
                    severity_V_2 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 0, waning_stage == "3-6 months"))$avg_ve_pred))/beta_coeffs$beta_V_2, 
                    severity_V_3 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 0, waning_stage == "6-12 months"))$avg_ve_pred))/beta_coeffs$beta_V_3,
                    severity_V_4 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 0, waning_stage == "12+ months"))$avg_ve_pred))/beta_coeffs$beta_V_4,
                    severity_1 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 1, waning_stage == "0-3 months"))$avg_ve_pred))/beta_coeffs$beta_1,
                    severity_2 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 1, waning_stage == "3-6 months"))$avg_ve_pred))/beta_coeffs$beta_2,
                    severity_3 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 1, waning_stage == "6-12 months"))$avg_ve_pred))/beta_coeffs$beta_3,
                    severity_4 = (1 - ((adj_waning_curves_severe %>% filter(prior_inf == 1, waning_stage == "12+ months"))$avg_ve_pred))/beta_coeffs$beta_4) 

write.csv(severity_coeffs, "data/clean-data/severity_coeffs.csv")

severity_coeffs_relabeled <- melt(severity_coeffs) %>% mutate(immunity_type = case_when(variable == "severity_0" ~ "immune-naive",
                                                                                                grepl("V", variable) ~ "vaccine-only",
                                                                                                TRUE ~ "hybrid/infection-only"),
                                                                      waning_group = case_when(grepl("1", variable) ~ "0-3 months",
                                                                                               grepl("2", variable) ~ "3-6 months",
                                                                                               grepl("3", variable) ~ "6-12 months",
                                                                                               grepl("4", variable) ~ "12+ months",
                                                                                               TRUE ~ NA)) %>% rename(severity_label = variable, severity_coeff = value)

write.csv(severity_coeffs_relabeled, "data/clean-data/severity_coeffs_relabeled.csv")
