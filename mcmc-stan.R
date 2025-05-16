######################################################
## PURPOSE: Script for fitting STAN models          ##
##          to RA serum sample data                 ##
######################################################

### Dependencies

library(rstan)
library(coda)

setwd("~/Documents/RA-Biomarker/")
source("~/Github/ra-biomarker/clean-data.R")
# source("~/Github/ra-biomarker/sens-clean-data.R") # sensitivity


## STAN Models

# Hyperparameters
a <- b <- rep(0, times = K)
R <- S <- diag(1e6, nrow = K, ncol = K) # beta/mu covariance hyperparameters

standata <- list(N = N, M = M, K = K, Y_obs = Y, L = L, U = U,
                 D_max = cens_max, D_obs = (1 - cens_max),
                 t = time, g = diagnosis, id = subj_id,
                 a = a, b = b, S = S, R = R)

# Truncation

stanmodel_trunc <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_trunc")
samples_trunc <- sampling(stanmodel_trunc, data = standata, iter = 5500, warmup = 500, 
                          chains = 4, thin = 5, check_data = FALSE, cores = 4)
mcmc_trunc <- do.call(cbind, rstan::extract(samples_trunc, 
                                            pars = c(paste0("gamma[", 1:K, "]"), 
                                                     paste0("kappa[", 1:K, "]")), 
                                            permuted = TRUE))

save(mcmc_trunc, file = "mcmc/mcmc_trunc.RData")
summary(mcmc_trunc)

# Censoring Above LOD

stanmodel_cens <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_cens")
samples_cens <- sampling(stanmodel_cens, data = standata, iter = 5500, warmup = 500, 
                         chains = 4, thin = 5, check_data = FALSE, cores = 4)
mcmc_cens <- do.call(cbind, rstan::extract(samples_cens, 
                                            pars = c(paste0("gamma[", 1:K, "]"), 
                                                     paste0("kappa[", 1:K, "]")), 
                                            permuted = TRUE))
save(mcmc_cens, file = "mcmc/mcmc_cens.RData")
summary(mcmc_cens)

# Binarize Detection

source("~/Github/ra-biomarker/sens-clean-data.R") # sensitivity

stanmodel_bin <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_bin")
samples_bin <- sampling(stanmodel_bin, data = standata, iter = 5500, warmup = 500, 
                          chains = 4, thin = 5, check_data = FALSE, cores = 4)
mcmc_bin <- do.call(cbind, rstan::extract(samples_bin, 
                                            pars = c(paste0("gamma[", 1:K, "]"), 
                                                     paste0("kappa[", 1:K, "]")), 
                                            permuted = TRUE))
save(mcmc_bin, file = "mcmc/mcmc_bin.RData")
summary(mcmc_bin)

### New Data

source("~/Github/ra-biomarker/clean-data-new.R")


a <- b <- rep(0, times = K)
R <- S <- diag(1e10, nrow = K, ncol = K) # beta/mu covariance hyperparameters

standata_new <- list(N = N, M = M, K = K, Y_obs = Y, L = L, U = U,
                     D_max = cens_max, D_obs = (1 - cens_max),
                     t = time, g = diagnosis, id = subj_id, 
                     a = a, b = b, S = S, R = R)


# Truncation

stanmodel_new_trunc <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_new_trunc")
samples_new_trunc <- sampling(stanmodel_new_trunc, data = standata_new, iter = 5500, warmup = 500,
                              chains = 4, thin = 5, check_data = FALSE, cores = 4)
mcmc_new_trunc <- do.call(cbind, rstan::extract(samples_trunc, 
                                                pars = c(paste0("gamma[", 1:K, "]"), 
                                                         paste0("kappa[", 1:K, "]")), 
                                                permuted = TRUE))
save(mcmc_new_trunc, file = "mcmc/mcmc_new_trunc.RData")
summary(mcmc_new_trunc)

# Censoring above LoD

stanmodel_new_cens <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_new_cens")
samples_new_cens <- sampling(stanmodel_new_cens, data = standata_new, iter = 5500, warmup = 500, 
                             chains = 4, thin = 5, check_data = TRUE, cores = 4)
mcmc_new_cens <- do.call(cbind, rstan::extract(samples_new_cens, 
                                              pars = c(paste0("gamma[", 1:K, "]"), 
                                                       paste0("kappa[", 1:K, "]")), 
                                              permuted = TRUE))
save(mcmc_new_cens, file = "mcmc/mcmc_new_cens.RData")
summary(mcmc_new_cens)

# Binarize Detection

source("~/Github/ra-biomarker/sens-clean-data-new.R")

stanmodel_new_bin <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_new_bin")
samples_new_bin <- sampling(stanmodel_new_bin, data = standata_new, iter = 5500, warmup = 500, 
                        chains = 4, thin = 5, check_data = TRUE, cores = 4)
mcmc_new_bin <- do.call(cbind, rstan::extract(samples_new_bin, 
                                              pars = c(paste0("gamma[", 1:K, "]"), 
                                                       paste0("kappa[", 1:K, "]")), 
                                              permuted = TRUE))
save(mcmc_new_bin, file = "mcmc/mcmc_new_bin.RData")
summary(mcmc_new_bin)

