######################################################
## PURPOSE: Script for fitting JAGS and STAN models ##
##          to RA serum sample data                 ##
## BY:      Kevin Josey                             ##
######################################################

### Dependencies

library(rjags)
load.module("glm")

source("~/Dropbox (Personal)/Projects/RA-Biomarker/code/clean-dat.R")
setwd("~/Dropbox (Personal)/Projects/RA-Biomarker/")

## JAGS Multivariate Model

# Hyperparameters
a <- b <- rep(0, times = K)
R <- S <- diag(1e-10, nrow = K, ncol = K) # beta/mu covariance hyperparameters
u <- 2
v <- rep(1e5, K)

X <- cbind(fem, bage, famhx)

# Correlated Intercepts

jagsDat <- list(N = N, M = M, K = K, Y = logY, t = t, id = study_id,
                a = a, b = b, S = S, R = R, u = u, v = v,
                race = nw, sex = fem, age = bage)

mod <- jags.model(file = "code/jags/change-point-a.jags", data = jagsDat, n.adapt = 100000, n.chains = 1)
mcmc <- coda.samples(mod, variable.names = c("beta", "kappa", "mu", "tau_e"), n.iter = 100000, thin = 100, na.rm = TRUE)

save(mcmc, file = "mcmc/mcmc_a.RData")
