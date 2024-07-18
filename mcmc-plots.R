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
mcmc_results <- read.csv("simulation-results/mcmc-chain-100sims-no-priors.csv")[,-1]
observed_data <- read.csv("data/clean-data/severe_cases_inc.csv", check.names = FALSE)[,-1] 

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

#Check calibration
final_params <- c(0.0001,0.0005, 0.005,0.05,0.05,1,1,1.2)

#final_params <- as.numeric(as.vector(mcmc_results[39,1:8]))
severity_param <- as.numeric(final_params[1:5])
initial_conditions <- model_pop_init(severity_param)
contact_matrix_adj <- contact_matrix_init(initial_conditions)
sim_with_final_params <- prediction(initial_conditions, final_params, contact_matrix_adj)


plot (R1 + R2 + R3 + R4 + R5 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "R")# ylim = c(0, 2000000)) 

plot (S_11 + S_12 + S_13 + S_14 + S_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S1", ylim = c(0, 2000000)) 
plot (E_11 + E_12 + E_13 + E_14 + E_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E1", ylim = c(0, 400000)) 
plot (I_11 + I_12 + I_13 + I_14 + I_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I1", ylim = c(0, 500000)) 
plot (H_11 + H_12 + H_13 + H_14 + H_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H1", ylim = c(0, 60)) 

plot (S_21 + S_22 + S_23 + S_24 + S_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S2") #, ylim = c(0, 2000000)) 
plot (E_21 + E_22 + E_23 + E_24 + E_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E2") #, ylim = c(0, 400000)) 
plot (I_21 + I_22 + I_23 + I_24 + I_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I2") #, ylim = c(0, 500000)) 
plot (H_21 + H_22 + H_23 + H_24 + H_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H2") #, ylim = c(0, 60)) 

plot (S_31 + S_32 + S_33 + S_34 + S_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S3", ylim = c(0, 2000000)) 
plot (E_31 + E_32 + E_33 + E_34 + E_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E3", ylim = c(0, 400000)) 
plot (I_31 + I_32 + I_33 + I_34 + I_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I3", ylim = c(0, 500000)) 
plot (H_31 + H_32 + H_33 + H_34 + H_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H3", ylim = c(0, 60)) 

plot (S_41 + S_42 + S_43 + S_44 + S_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S4", ylim = c(0, 2000000)) 
plot (E_41 + E_42 + E_43 + E_44 + E_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E4", ylim = c(0, 400000)) 
plot (I_41 + I_42 + I_43 + I_44 + I_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I4", ylim = c(0, 500000)) 
plot (H_41 + H_42 + H_43 + H_44 + H_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H4", ylim = c(0, 60)) 

plot (S_V_11 + S_V_12 + S_V_13 + S_V_14 + S_V_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S_V_1", ylim = c(0, 2000000)) 
plot (E_V_11 + E_V_12 + E_V_13 + E_V_14 + E_V_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E_V_1", ylim = c(0, 400000)) 
plot (I_V_11 + I_V_12 + I_V_13 + I_V_14 + I_V_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I_V_1", ylim = c(0, 500000)) 
plot (H_V_11 + H_V_12 + H_V_13 + H_V_14 + H_V_15 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H_V_1", ylim = c(0, 60)) 

plot (S_V_21 + S_V_22 + S_V_23 + S_V_24 + S_V_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S_V_2", ylim = c(0, 2000000)) 
plot (E_V_21 + E_V_22 + E_V_23 + E_V_24 + E_V_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E_V_2", ylim = c(0, 400000)) 
plot (I_V_21 + I_V_22 + I_V_23 + I_V_24 + I_V_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I_V_2", ylim = c(0, 500000)) 
plot (H_V_21 + H_V_22 + H_V_23 + H_V_24 + H_V_25 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H_V_2", ylim = c(0, 60)) 

plot (S_V_31 + S_V_32 + S_V_33 + S_V_34 + S_V_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S_V_3", ylim = c(0, 2000000)) 
plot (E_V_31 + E_V_32 + E_V_33 + E_V_34 + E_V_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E_V_3", ylim = c(0, 400000)) 
plot (I_V_31 + I_V_32 + I_V_33 + I_V_34 + I_V_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I_V_3", ylim = c(0, 500000)) 
plot (H_V_31 + H_V_32 + H_V_33 + H_V_34 + H_V_35 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H_V_3", ylim = c(0, 60)) 

plot (S_V_41 + S_V_42 + S_V_43 + S_V_44 + S_V_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "S_V_4", ylim = c(0, 2000000)) 
plot (E_V_41 + E_V_42 + E_V_43 + E_V_44 + E_V_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "E_V_4", ylim = c(0, 400000)) 
plot (I_V_41 + I_V_42 + I_V_43 + I_V_44 + I_V_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "I_V_4", ylim = c(0, 500000)) 
plot (H_V_41 + H_V_42 + H_V_43 + H_V_44 + H_V_45 ~ time, data = as.data.frame(sim_with_final_params), type='b', col = 'blue', xlab = "Time (days)", ylab = "H_V_4", ylim = c(0, 60)) 

-log(score(sim_with_final_params, observed_data))

#Plot
plot_data <- melt(sim_with_final_params, id = "week")
ggplot(data = plot_data, aes(x = week, y = value, color = variable)) +
  geom_line() +
  ylim(0, 500) +
  ylab("Severe Case Incidence (per 100,000)")+
  xlab("Time")+
  ggtitle("Calibrated Severe Case Incidence Fit by age (low SSE)")

ggplot(data = melt(observed_data, id = "week"), aes(x = as.Date(week), y = value, color = variable)) +
  geom_line() +
  ylim(0, 500) +
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

sse <- data.frame(sse = exp(-mcmc_results$V9),
                  sims = c(1:101))

ggplot(data = sse, aes(sims, sse)) + 
  geom_line() + geom_point()+
  xlab("Number of Sims") +
  ylab("SSE") +
  ggtitle("SSE scores over mcmc calibration")

count_same <- 0
for (i in c(2: nrow(mcmc_results))){
  if(setequal(mcmc_results[i,], mcmc_results[i-1,])){
    count_same <- count_same + 1
    print(i)
  }
  
}

acceptance_prob <- c()
for (i in c(2: nrow(mcmc_results))){
  acceptance_prob <- append(acceptance_prob, min(exp(mcmc_results[i,9] - mcmc_results[i-1,9]), 1) - 0.7)
}

min(acceptance_prob)
max(acceptance_prob)

