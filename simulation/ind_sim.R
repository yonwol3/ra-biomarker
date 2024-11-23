################################################
##### Simulation study for one outcome ##########
###############################################
library(tidyverse)
library(rstan)

# Ei ~ N(0,sigma_e^2)
# alphai ~ N(0,sigma_b^2)
# timeij ~ unif (a,b)
# We have 20 subject IDs and 1:1 ratio for treatment allocation



# Parameters
sigma_b <- 1  # Standard deviation for random intercepts
sigma_e <- 2  # Standard deviation for error terms
n <- 200      # Total number of observations
m <- 50       # Number of subjects
c <- n/m     # Observations per subject
kappa <- -8   # Change point parameter
a <- -20      # Uniform distribution lower bound
b <- 5        # Uniform distribution upper bound
b_1 <- 0.1    # Time parameter
b_2 <- 0.6    # Treatment parameter
b_3 <- 0.3    # Interaction term

true_beta3 <- 0.3
true_kappa <- -8
n_sim <- 1000



set.seed(123)  # For reproducibility



# Create empty data frames to store results
results <- data.frame(
  sim_id = 1:n_sim,
  beta3_mean = NA, beta3_median = NA, beta3_lower = NA, beta3_upper = NA,
  kappa_mean = NA, kappa_median = NA, kappa_lower = NA, kappa_upper = NA
)

# Power, bias, MSE, coverage placeholders for each parameter
metrics <- data.frame(
  beta3_bias = NA, beta3_mse = NA, beta3_power = NA, beta3_coverage = NA,
  kappa_bias = NA, kappa_mse = NA, kappa_power = NA, kappa_coverage = NA
)

# Simulate 10,000 times
for (i in 1:n_sim) {
  
  id <- rep(1:m, each=c)
  time <- runif(n, a, b)
  ccp3 <- as.numeric(id <= 10)  # 1 for first 10 subjects, 0 for others
  alpha <- rep(rnorm(m, 0, sigma_b), each=c)  # Random intercepts
  
  # Outcome calculation
  error_term <- rnorm(n, 0, sigma_e)
  outcome <- alpha + b_1 * time + b_2 * ccp3 + ifelse(time > kappa & ccp3 == 1, b_3 * (time - kappa), 0) + error_term
  
  # Data frame
  df <- data.frame(observation = 1:n, id = id, time = time, ccp3 = ccp3, alpha = alpha, outcome = outcome)
  
  
  # Fit the Stan model
   stan_dat<- list(N=nrow(df),
                   M= length(unique(df$id)),
                   Y=df$outcome, # outcome assumed to be in log scale 
                   t=df$time,
                   g=df$ccp3,
                   id=df$id,
                   a=0, 
                   r= 10^6)
  
  stan_fit<- stan(data = stan_dat, file = "stan/ind_sim.stan", 
                  chains = 1, iter = 4000, warmup = 2000)
  
  # Extract summary statistics
  fit_summary <- summary(stan_fit, pars = c("beta3", "kappa"))$summary
  
  # Save the mean, median, and credible intervals for beta3 and kappa
  results$beta3_mean[i] <- fit_summary["beta3", "mean"]
  results$beta3_median[i] <- fit_summary["beta3", "50%"]
  results$beta3_lower[i] <- fit_summary["beta3", "2.5%"]
  results$beta3_upper[i] <- fit_summary["beta3", "97.5%"]
  
  results$kappa_mean[i] <- fit_summary["kappa", "mean"]
  results$kappa_median[i] <- fit_summary["kappa", "50%"]
  results$kappa_lower[i] <- fit_summary["kappa", "2.5%"]
  results$kappa_upper[i] <- fit_summary["kappa", "97.5%"]
}

# Calculate metrics (Power, Bias, MSE, and Coverage Probability)
metrics <- results %>%
  mutate(
    # Bias: Mean - True value
    beta3_bias = beta3_mean - true_beta3,
    kappa_bias = kappa_mean - true_kappa,
    
    # MSE: (Mean - True value)^2
    beta3_mse = (beta3_mean - true_beta3)^2,
    kappa_mse = (kappa_mean - true_kappa)^2,
    
    # Power: Proportion of times the 95% credible interval excludes 0
    beta3_power = ifelse(beta3_lower > 0 | beta3_upper < 0, 1, 0),
    kappa_power = ifelse(kappa_lower > 0 | kappa_upper < 0, 1, 0),
    
    # Coverage Probability: Proportion of times the 95% credible interval contains the true value
    beta3_coverage = ifelse(beta3_lower <= true_beta3 & beta3_upper >= true_beta3, 1, 0),
    kappa_coverage = ifelse(kappa_lower <= true_kappa & kappa_upper >= true_kappa, 1, 0)
  )

# Summarize results (mean across simulations for each metric)
final_results <- data.frame(
  beta3_bias = mean(metrics$beta3_bias),
  beta3_mse = mean(metrics$beta3_mse),
  beta3_power = mean(metrics$beta3_power),
  beta3_coverage = mean(metrics$beta3_coverage),
  
  kappa_bias = mean(metrics$kappa_bias),
  kappa_mse = mean(metrics$kappa_mse),
  kappa_power = mean(metrics$kappa_power),
  kappa_coverage = mean(metrics$kappa_coverage)
)

# Display the final summarized results
print(final_results)


