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
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)

standata <- list(N = N, M = M, K = K, Y_obs = logY, L = L, U = U,
                 D_min = cens_min, D_max = cens_max, D_obs = (1 - cens_max)*(1 - cens_min),
                 t = time, g = diagnosis, id = subj_id, a = a, b = b, S = S, R = R)

stanmodel_a <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-a.stan", model_name = "stanmodel_a")
samples_a <- sampling(stanmodel_a, data = standata, iter = 20000, warmup = 5000, chains = 1, thin = 15, check_bata = FALSE)
mcmc_a <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_a, pars = c(paste0("gamma[", 1:K, "]"),
                                                                            paste0("kappa[", 1:K, "]")),
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_a, file = "mcmc/mcmc_a.RData")
summary(mcmc_a)

stanmodel_b <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-b.stan", model_name = "stanmodel_b")
samples_b <- sampling(stanmodel_b, data = standata, iter = 20000, warmup = 5000, chains = 1, thin = 15, check_bata = FALSE)
mcmc_b <- do.call(mcmc.list, plyr::alply(rstan::extract(samples_b, pars = c(paste0("gamma[", 1:K, "]"), 
                                                                            paste0("kappa[", 1:K, "]")), 
                                                        permuted = FALSE), 2, coda::mcmc))
save(mcmc_b, file = "mcmc/mcmc_b.RData")
summary(mcmc_b)
