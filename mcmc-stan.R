######################################################
## PURPOSE: Script for fitting STAN models          ##
##          to RA serum sample data                 ##
######################################################

### Dependencies

library(rstan)
library(coda)

source("~/Github/ra-biomarker/clean-data.R")
setwd("~/Documents/RA-Biomarker/")

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
samples_trunc <- sampling(stanmodel_trunc, data = standata, iter = 30000, warmup = 5000, 
                          chains = 1, thin = 25, check_data = FALSE)
mcmc_trunc <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_trunc, 
                                                            pars = c(paste0("gamma[", 1:K, "]"), 
                                                                     paste0("kappa[", 1:K, "]")), 
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_trunc, file = "mcmc/mcmc_trunc.RData")
summary(mcmc_trunc)

# Censoring Above LOD

stanmodel_cens <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_cens")
samples_cens <- sampling(stanmodel_cens, data = standata, iter = 30000, warmup = 5000, 
                         chains = 1, thin = 25, check_data = FALSE)
mcmc_cens <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_cens, 
                                                           pars = c(paste0("gamma[", 1:K, "]"), 
                                                                    paste0("kappa[", 1:K, "]")), 
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_cens, file = "mcmc/mcmc_cens.RData")
summary(mcmc_cens)

### New Data

source("~/Github/ra-biomarker/clean-data-new.R")

for (bio_dat in bio_couple_data) {
  
  a <- b <- rep(0, times = bio_dat$K)
  R <- S <- diag(1e6, nrow = bio_dat$K, ncol = bio_dat$K) # beta/mu covariance hyperparameters
  
  standata_new <- with(bio_dat,
                       list(N = N, M = M, K = K, Y_obs = Y, L = L, U = U,
                            D_max = cens_max, D_obs = (1 - cens_max),
                            t = time, g = diagnosis, id = subj_id, 
                            a = a, b = b, S = S, R = R))
  
  stanmodel_new_trunc <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", 
                                    model_name = paste0("stanmodel_new_trunc_", bio_dat$biomarker))
  samples_new_trunc <- sampling(stanmodel_new_trunc, data = standata_new, iter = 30000, warmup = 5000, 
                                chains = 1, thin = 25, check_data = FALSE)
  mcmc_new_trunc <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_new_trunc,
                                                                  pars = c(paste0("gamma[", 1:bio_dat$K, "]"), 
                                                                           paste0("kappa[", 1:bio_dat$K, "]")), 
                                                                  permuted = FALSE), 2, coda::mcmc))
  save(mcmc_new_trunc, file = paste0("mcmc/mcmc_new_trunc_", bio_dat$biomarker, ".RData"))
  summary(mcmc_new_trunc)
  
  
}

