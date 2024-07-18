########################################################################################################################
#Title: MCMC Calibration
#Author: Hailey Park
#Date: April 16th, 2024
########################################################################################################################

prediction <- function (initial_conditions, params, contact_matrix_adj) {

  #Incorporate the contact matrix factors into the set of params
  add_params <- c(params, contact_matrix_adj)

  #Run model
  results <- lsoda(initial_conditions, time, seir_model, add_params)
  
  #Reformat results into weekly incidence estimates by age group (separate columns)
  reformatted_results <- as.data.frame(results) %>% select(time, new_severe_cases1:new_severe_cases5) %>%
    mutate(week = time %/% 7) %>% group_by(week) %>% summarise(across(c(new_severe_cases1:new_severe_cases5), sum)) %>%
    mutate(age_0_17_inc = new_severe_cases1/total_pop_by_age$total_pop[1] * 100000,
           age_18_49_inc = new_severe_cases2/total_pop_by_age$total_pop[2] * 100000,
           age_50_64_inc = new_severe_cases3/total_pop_by_age$total_pop[3] * 100000,
           age_65_74_inc = new_severe_cases4/total_pop_by_age$total_pop[4] * 100000,
           age_75plus_inc = new_severe_cases5/total_pop_by_age$total_pop[5] * 100000) %>%
    select(-c(new_severe_cases1:new_severe_cases5))

  return(reformatted_results)
}


# The score is the L² distance of the simulated inc predictions from the observed data.
score <- function (sim_pred, data) {

  age_0_17 <- sum(data$`0-17 years` - sim_pred$age_0_17_inc)^2
  age_18_49 <- sum(data$`18-49 years` - sim_pred$age_18_49_inc)^2
  age_50_64 <- sum(data$`50-64 years` - sim_pred$age_50_64_inc)^2
  age_65_74 <- sum(data$`65-74 years` - sim_pred$age_65_74_inc)^2
  age_75_plus <- sum(data$`75+ years` - sim_pred$age_75plus_inc)^2
  
  total <- sum(age_0_17, age_18_49, age_50_64, age_65_74, age_75_plus)
  return(total)
}


# The likelihood is the negative log score of the average simulated incidence.
likelihood <- function (params) {
  
  #Parameters
  severity_by_age <- c(params[1], params[2], params[3], params[4], params[5])
  nonsevere_waning_rate_vaccine <- params[6]
  nonsevere_waning_rate_hybrid <- params[7]
  lambda <- params[8]
  
  #Model initialization
  initial_conditions <- model_pop_init(severity_by_age)
  
  #If initial conditions have negative totals in the susceptible compartments (if severity parameters are too small),
  # return null
  if(any(is.na(initial_conditions))) {
    return(-Inf)
  }
  
  #Calculate contact matrix adjustment factors
  contact_matrix_adj <- contact_matrix_init(initial_conditions)
  
  #Get predictions
  sim_pred <- prediction(initial_conditions, params, contact_matrix_adj)
  
  
  return(-log(score(sim_pred, observed_data)))
  
}

# Prior distribution
prior <- function(par_initial){
  
  severity_0_17 <- par_initial[1]
  severity_18_49 <- par_initial[2]
  severity_50_64 <- par_initial[3]
  severity_65_74 <- par_initial[4]
  severity_75plus <- par_initial[5]
  pe_nonsevere_rate_vacc = par_initial[6]
  pe_nonsevere_rate_hybrid = par_initial[7]
  lambda = par_initial[8]
  
  severity_0_17_prior <- dunif(severity_0_17, min=0.000005, max=0.0005, log=TRUE)
  severity_18_49_prior <- dunif(severity_18_49, min=0.0001, max=0.001, log=TRUE)
  severity_50_64_prior <- dunif(severity_50_64, min=0.001, max=0.01, log=TRUE)
  severity_65_74_prior <- dunif(severity_65_74, min=0.005, max=0.1, log=TRUE)
  severity_75plus_prior <- dunif(severity_75plus, min=0.005, max=0.1, log=TRUE)
  pe_nonsevere_rate_vacc_prior <- dlnorm(pe_nonsevere_rate_vacc, meanlog=0.1, sdlog=0.8, log=TRUE)
  pe_nonsevere_rate_hybrid_prior <- dlnorm(pe_nonsevere_rate_hybrid, meanlog=0.1, sdlog=0.8, log=TRUE)
  lambda_prior <- dunif(lambda, min=0.5, max=2, log=TRUE)
  
  return(-(severity_0_17_prior + severity_18_49_prior + severity_50_64_prior + severity_65_74_prior + severity_75plus_prior +
           pe_nonsevere_rate_vacc_prior+pe_nonsevere_rate_hybrid_prior+
           lambda_prior))
}

posterior <- function(param){
  likel <- likelihood(param)
  priors <- prior(param)
  print(paste0("Likelihood: ", likel))
  print(paste0("Priors: ", priors))
  return (likel + priors)
 # return (likelihood(param) + prior(param))
}


# Choosing a new parameter value by sampling from a multivariate normal distribution
# centered at the current value, that is called the proposal function. The SDs 
# roughly correspond to the step-size of the chain for each parameter.
proposalfunction <- function(param){
  
  severity_0_17 <- min(max(rnorm(1, mean=param[1], sd=0.000005), 0.000005), 0.0005)
  severity_18_49 <- min(max(rnorm(1, mean=param[2], sd=0.00005),0.0001), 0.001)
  severity_50_64 <- min(max(rnorm(1, mean=param[3], sd=0.0005), 0.001), 0.01)
  severity_65_74 <- min(max(rnorm(1, mean=param[4], sd=0.005), 0.005), 0.1)
  severity_75plus <- min(max(rnorm(1, mean=param[5], sd=0.005), 0.005), 0.1)
  pe_nonsevere_rate_vacc = min(max(rnorm(1, mean=param[6], sd=0.1), 0.3), 4)
  pe_nonsevere_rate_hybrid = min(max(rnorm(1, mean=param[7], sd=0.1), 0.3), 4)
  lambda = min(max(rnorm(1, mean=param[8], sd=0.01), 0.5), 2)
  
  return(c(severity_0_17, severity_18_49, severity_50_64, severity_65_74, severity_75plus,
           pe_nonsevere_rate_vacc, pe_nonsevere_rate_hybrid, lambda))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,8))
  chain[1,] = startvalue
  print(startvalue)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Proposal: ")
    print(proposal)
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    random_probab <- runif(1)
    
    # print(paste0("Probability: ", probab))
    # print(paste0("Random prob: ", random_probab))
    
    if (random_probab < probab){
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
    
    print(paste0("Iteration ", i, " Chain: "))
    print(chain[i+1,])
  }
  return(chain)
}
