######################################################
## PURPOSE: Script for fitting STAN models          ##
##          to RA serum sample data                 ##
######################################################

### Dependencies

library(rstan)
library(coda)

source("~/Github/ra-biomarker/clean-data.R")
setwd("~/Dropbox/Projects/RA-Biomarker/")

## STAN Models

# Hyperparameters
a <- b <- rep(0, times = K)
R <- S <- diag(1e6, nrow = K, ncol = K) # beta/mu covariance hyperparameters

standata <- list(N = N, M = M, K = K, Y_obs = logY, 
                 L = minY, U = maxY,
                 D_min = cens_min, D_max = cens_max, 
                 D_obs = (1 - cens_max)*(1 - cens_min),
                 t = time, g = diagnosis, id = subj_id,
                 a = a, b = b, S = S, R = R)

stanmodel_c <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-a.stan", model_name = "stanmodel_a")
samples_c <- sampling(stanmodel_c, data = standata, iter = 20000, warmup = 5000, chains = 1, thin = 15, check_data = FALSE)
mcmc_c <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_c, pars = c(paste0("gamma[", 1:K, "]"),
                                                                            paste0("kappa[", 1:K, "]")),
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_c, file = "mcmc/mcmc_a.RData")
summary(mcmc_c)

stanmodel_d <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-b.stan", model_name = "stanmodel_b")
samples_d <- sampling(stanmodel_d, data = standata, iter = 20000, warmup = 5000, chains = 1, thin = 15, check_data = FALSE)
mcmc_c <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_c, pars = c(paste0("gamma[", 1:K, "]"), 
                                                                            paste0("kappa[", 1:K, "]")), 
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_d, file = "mcmc/mcmc_b.RData")
summary(mcmc_d)
