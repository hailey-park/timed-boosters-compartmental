########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################

rm(list=ls())
gc()

#Load libraries
library(dplyr)
library(ggplot2)
library(tibble)
library(lubridate)
library(here)
library(data.table)
library(tidyverse)
library(deSolve) 
library(reshape2)

options(dplyr.summarise.inform = FALSE)

# Observed data
observed_data <- read.csv("data/clean-data/severe_cases_inc.csv", check.names = FALSE)[,-1] 
  
#Set initial params
starting_params <- c(severity_0_17 = 0.0001,
                     severity_18_49 = 0.0005,
                     severity_50_64 = 0.005,
                     severity_65_74 = 0.05,
                     severity_75plus = 0.05,
                     nonsevere_waning_rate_vaccine = 1,
                     nonsevere_waning_rate_hybrid = 1,
                     lambda = 1)

# Reading model scripts
source(here::here("seir-model.R"))
source(here::here("model-pop-init.R"))
source(here::here("mcmc-calibration-functions-no-priors.R"))

#Simulation vars
time <- seq(0, 494, by=1) # Number of days. (Jan 1, 2022 - May 10, 2023)
latent_period <- 3        # Average latent period length.
infectious_period <- 5    # Average infectious period length.
gamma <- 1/infectious_period       
delta <- 1/latent_period
N <- 39512223             # Total population size


#Set simulation size
sim_size <- 10000

#Run MCMC
chain = run_metropolis_MCMC(starting_params, sim_size)

# burnIn = 1000
# acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

#Set up folder structure to save simulation results
dir.create("simulation-results")

#Save results
write.csv(chain, "simulation-results/mcmc-chain-10000sims-no-priors-no-adj.csv")

# params["beta"] <- chain[sim_size, 2]
# params["gamma"] <- chain[sim_size, 1]
# times <- c(data$week)
# model.pred <- prediction(params,times)*params["N"]
# plot(incidence~week,data=data,type='b',col='red') 
# lines(times,model.pred,type='l')



