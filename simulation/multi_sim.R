##########################################
## multivariate simulation with stan
## change_point a (correlation on error)
## By: Yonatan Wolde
##########################################

library(tidyverse)
library(MASS)
library(rstan)
library(clusterGeneration)  # For generating positive definite covariance matrices

set.seed(2222)

#-------------------------#
# Parameter Specification
#-------------------------#

a <- -20      # Uniform distribution lower bound
b <- 5        # Uniform distribution upper bound
n <- 500      # Total number of observations
m <- 50       # Number of subjects
c <- n / m    # Observations per subject
k <- 6        # Number of biomarkers
sigma_0 <- 1  # Standard deviation for random intercepts
sigma_e <- 2  # Standard deviation for error terms

# Variance-covariance matrix for random intercepts (assuming uncorrelated)
var_cov_sigma_0 <- diag(sigma_0^2, k)  # Variance-covariance matrix for alpha

# Generate a random correlation matrix for errors
rho_e <- genPositiveDefMat(k, covMethod = "onion")$Sigma

# Create the covariance matrix for errors
var_cov_sigma_e <- diag(sigma_e, k) %*% rho_e %*% diag(sigma_e, k)

beta1 <- rnorm(k)         # Time parameters
beta2 <- rnorm(k)         # Treatment parameters
gamma <- rnorm(k)         # Interaction parameters
kappa <- runif(k, min = a, max = b)  # Change point parameters
mu <- rnorm(k)

#---------------------------#
# Simulation Study
#--------------------------# 

Y <- matrix(NA, nrow = n, ncol = k)

id <- rep(1:m, each = c)
time <- runif(n, a, b)
ccp3 <- as.numeric(id <= 20)
alpha <- mvrnorm(m, mu = mu, Sigma = var_cov_sigma_0)     # Simulate random intercepts
error <- mvrnorm(n, mu = rep(0, k), Sigma = var_cov_sigma_e)  # Simulate random errors

for (j in 1:k) {
  Y[, j] <- rep(alpha[, j], each = c) + beta1[j] * time + beta2[j] * ccp3 +
    ifelse(time > kappa[j] & ccp3 == 1, gamma[j] * (time - kappa[j]), 0) + error[, j]
}

dat <- as.data.frame(cbind(Y, id = id, time = time, ccp3 = ccp3))

stan_dat <- list(
  N = nrow(dat),
  M = length(unique(dat$id)),
  K = k,
  id = dat$id,
  t = dat$time,
  g = dat$ccp3,
  a = rep(0, k),
  b = rep(0, k),
  R = diag(1e6, nrow = k),
  S = diag(1e6, nrow = k),
  Y = as.matrix(dat[, 1:k])
)

stan_fit <- stan(
  data = stan_dat,
  file = "stan/multi_sim.stan",
  chains = 1,
  iter = 4000,
  warmup = 2000
)

pars=c(paste0("gamma","[",1:k,"]"),paste0("kappa","[",1:k,"]")) # specific parameters for our research question

summary_stanfit<-summary(stan_fit, pars=pars)$summary
summary_stanfit
rstan::traceplot(stan_fit, pars = pars)
