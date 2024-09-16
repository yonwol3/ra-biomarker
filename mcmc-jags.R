######################################################
## PURPOSE: Script for fitting JAGS and STAN models ##
##          to RA serum sample data                 ##
## BY:      Kevin Josey                             ##
######################################################

### Dependencies

library(rjags)
library(coda)
load.module("glm")

source("~/Github/ra-biomarker/clean-dat.R")
setwd("~/Dropbox/Projects/RA-Biomarker/")

## JAGS Multivariate Model

# Hyperparameters
a <- b <- rep(0, times = K)
R <- S <- diag(1e-6, nrow = K, ncol = K) # beta/mu covariance hyperparameters
u <- 2
v <- rep(1e-6, K)

jagsDat <- list(N = N, M = M, K = K, Y = logY, t = t, id = study_id,
                a = a, b = b, S = S, R = R, u = u, v = v,
                race = nw, sex = fem, age = bage)

# Uncorrelated Intercepts

mod_a <- jags.model(file = "~/Github/ra-biomarker/jags/change-point-a.jags", data = jagsDat, n.adapt = 100000, n.chains = 1)
mcmc_a <- coda.samples(mod_a, variable.names = c("beta", "kappa", "mu", "tau_e"), n.iter = 100000, thin = 100, na.rm = TRUE)

save(mcmc_a, file = "mcmc/mcmc_a.RData") 

# Correlated Intercepts

mod_b <- jags.model(file = "~/Github/ra-biomarker/jags/change-point-b.jags", data = jagsDat, n.adapt = 100000, n.chains = 1)
mcmc_b <- coda.samples(mod_b, variable.names = c("beta", "kappa", "mu", "tau_e"), n.iter = 100000, thin = 100, na.rm = TRUE)

save(mcmc_b, file = "mcmc/mcmc_b.RData")