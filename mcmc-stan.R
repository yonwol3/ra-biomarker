######################################################
## PURPOSE: Script for fitting JAGS and STAN models ##
##          to RA serum sample data                 ##
## BY:      Kevin Josey                             ##
######################################################

### Dependencies

library(rstan)
library(coda)

source("~/Github/ra-biomarker/clean-dat.R")
setwd("~/Dropbox/Projects/RA-Biomarker/")

## STAN Models

# Hyperparameters
a <- b <- rep(0, times = K)
R <- S <- diag(1e6, nrow = K, ncol = K) # beta/mu covariance hyperparameters

standata <- list(N = N, M = M, K = K, Y = logY, 
                 t = t, g = diagnosis, id = subj_id,
                 a = a, b = b, S = S, R = R)

stanmodel_c <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-c.stan", model_name="stanmodel_c")
samples_c <- sampling(stanmodel_c, data = standata, iter = 20000, warmup = 5000, chains = 1, check_data = FALSE)
mcmc_c <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_c, pars = c(paste0("gamma[", 1:K, "]"),
                                                                            paste0("kappa[", 1:K, "]")), 
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_c, file = "mcmc/mcmc_c.RData")
summary(mcmc_c)

stanmodel_d <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-d.stan", model_name="stanmodel_d")
samples_d <- sampling(stanmodel_d, data = standata, iter = 20000, warmup = 5000, chains = 1, check_data = FALSE)
mcmc_c <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_c, pars = c(paste0("gamma[", 1:K, "]"), 
                                                                            paste0("kappa[", 1:K, "]")), 
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_d, file = "mcmc/mcmc_d.RData")
summary(mcmc_d)
