######################################################
## PURPOSE: Script for fitting STAN models          ##
##          to RA serum sample data                 ##
######################################################

### Dependencies

library(rstan)
library(coda)

setwd("~/Documents/RA-Biomarker/")
source("~/Github/ra-biomarker/clean-data-A.R")

## STAN Models

# Hyperparameters
a <- b <- rep(0, times = K)
R <- S <- diag(1e8, nrow = K, ncol = K) # beta/mu covariance hyperparameters

standata_A <- list(N = N, M = M, K = K, Y = Y, L = L, U = U, D = D,
                 t = time, g = diagnosis, id = study_id, a = a, b = b, S = S, R = R)

# Truncation at LoD
stanmodel_trunc_A <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_trunc_A")
samples_trunc_A <- sampling(stanmodel_trunc_A, data = standata_A, iter = 20000, warmup = 10000, 
                          chains = 5, thin = 10, check_data = FALSE, cores = 5,
                          control = list(max_treedepth = 15))
mcmc_trunc_A <- do.call(cbind, rstan::extract(samples_trunc_A, 
                                              pars = c(paste0("gamma[", 1:K, "]"), 
                                                       paste0("delta[", 1:K, "]")), 
                                              permuted = TRUE))
save(mcmc_trunc_A, file = "mcmc/mcmc_trunc_A.RData")
summary(mcmc_trunc_A)

# Censoring Above LoD
stanmodel_cens_A <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_cens_A")
samples_cens_A <- sampling(stanmodel_cens_A, data = standata_A, iter = 20000, warmup = 10000, 
                         chains = 5, thin = 10, check_data = FALSE, cores = 5, 
                         control = list(max_treedepth = 15))
mcmc_cens_A <- do.call(cbind, rstan::extract(samples_cens_A, 
                                             pars = c(paste0("gamma[", 1:K, "]"), 
                                                      paste0("delta[", 1:K, "]")), 
                                            permuted = TRUE))
save(mcmc_cens_A, file = "mcmc/mcmc_cens_A.RData")
summary(mcmc_cens_A)

### Sample B

source("~/Github/ra-biomarker/clean-data-B.R")

a <- b <- rep(0, times = K)
R <- S <- diag(1e8, nrow = K, ncol = K) # beta/mu covariance hyperparameters

standata_B <- list(N = N, M = M, K = K, Y = Y, L = L, U = U, D = D,
                     t = time, g = diagnosis, id = study_id, 
                     a = a, b = b, S = S, R = R)

# Truncation at LoD
stanmodel_trunc_B <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_trunc_B")
samples_trunc_B <- sampling(stanmodel_trunc_B, data = standata_B, iter = 20000, warmup = 10000,
                            chains = 5, thin = 10, check_data = FALSE, cores = 5,
                            control = list(max_treedepth = 15))
mcmc_trunc_B <- do.call(cbind, rstan::extract(samples_trunc_B, 
                                              pars = c(paste0("gamma[", 1:K, "]"), 
                                                       paste0("delta[", 1:K, "]")), 
                                              permuted = TRUE))
save(mcmc_trunc_B, file = "mcmc/mcmc_trunc_B.RData")
summary(mcmc_trunc_B)

# Censoring above LoD
stanmodel_cens_B <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_cens_B")
samples_cens_B <- sampling(stanmodel_cens_B, data = standata_B, iter = 20000, warmup = 10000, 
                           chains = 5, thin = 10, check_data = TRUE, cores = 5,
                           control = list(max_treedepth = 15))
mcmc_cens_B <- do.call(cbind, rstan::extract(samples_cens_B, 
                                             pars = c(paste0("gamma[", 1:K, "]"), 
                                                      paste0("delta[", 1:K, "]")), 
                                             permuted = TRUE))
save(mcmc_cens_B, file = "mcmc/mcmc_cens_B.RData")
summary(mcmc_cens_B)

