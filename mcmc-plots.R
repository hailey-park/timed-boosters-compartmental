###################################################################################################
#Title: MCMC plots
#Author: Hailey Park
#Date: July 15, 2024
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
library(deSolve) 
library(reshape2)

#Read data
mcmc_results <- read.csv("simulation-results/mcmc-chain-10000sims-no priors.csv")[,-1]
observed_data <- read.csv("data/clean-data/severe_cases_inc.csv", check.names = FALSE)[,-1] 

# Reading model scripts
source(here::here("seir-model.R"))
source(here::here("model-pop-init.R"))
source(here::here("mcmc-calibration-functions.R"))

#Simulation vars
time <- seq(0, 494, by=1) # Number of days. (Jan 1, 2022 - May 10, 2023)
latent_period <- 3        # Average latent period length.
infectious_period <- 5    # Average infectious period length.
gamma <- 1/infectious_period       
delta <- 1/latent_period
N <- 39512223             # Total population size

#Check calibration
final_params <- as.numeric(as.vector(mcmc_results[nrow(mcmc_results),]))
severity_param <- as.numeric(final_params[1:5])
initial_conditions <- model_pop_init(severity_param)
contact_matrix_adj <- contact_matrix_init(initial_conditions)
sim_with_final_params <- prediction(initial_conditions, final_params, contact_matrix_adj)


#Plot
plot_data <- melt(sim_with_final_params, id = "week")
ggplot(data = plot_data, aes(x = week, y = value, color = variable)) +
  geom_line() +
  ylim(0, 400) +
  ylab("Severe Case Incidence (per 100,000)")+
  xlab("Time")+
  ggtitle("Calibrated Severe Case Incidence Fit by age")

ggplot(data = melt(observed_data, id = "week"), aes(x = as.Date(week), y = value, color = variable)) +
  geom_line() +
  ylim(0, 400) +
  ylab("Severe Case Incidence (per 100,000)")+
  xlab("Time")+
  ggtitle("Observed Severe Case Incidence by age")
########################################################


#Plot parameter values tested over MCMC runs

plot(mcmc_results$V1)
plot(mcmc_results$V2)
plot(mcmc_results$V3)
plot(mcmc_results$V4)
plot(mcmc_results$V5)
plot(mcmc_results$V6)
plot(mcmc_results$V7)
plot(mcmc_results$V8)

