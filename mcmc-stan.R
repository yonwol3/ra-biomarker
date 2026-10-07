######################################################
## PURPOSE: fit STAN models to RA serum data        ##
##                                                  ##
## PRIMARY:     change-point-cens-corr.stan         ##
##              censored (data augmentation) +      ##
##              LKJ correlation, on the log scale   ##
##              with weakly-informative slopes.     ##
## SENSITIVITY: change-point-cens-log.stan          ##
##              SAME log scale, SAME priors, SAME   ##
##              censoring; DIAGONAL residual cov.   ##
##              Differs from the primary ONLY in    ##
##              the correlation structure.          ##
######################################################

### Dependencies

library(rstan)
library(coda)

setwd("~/Documents/RA-Biomarker/")
source("~/Github/ra-biomarker/clean-data-A.R")

set.seed(42)

## STAN Models

### Sample A

# Hyperparameters
a <- b <- rep(0, times = K)
R    <- diag(1e6, nrow = K, ncol = K)   # theta (intercept) prior covariance
S    <- diag(1e6, nrow = K, ncol = K)   # beta/gamma prior covariance (sensitivity)
S_cc <- diag(4, nrow = K, ncol = K)   # weakly-informative slope prior (log scale) for the primary

# Data augmentation indices for the censored cells (D == 1)
cens <- which(D == 1, arr.ind = TRUE)
cens <- cens[order(cens[,1], cens[,2]),]

standata_cc_A <- list(N = N, M = M, K = K, Y = log(Y), U = log(U),
                      t = time, g = diagnosis, id = study_id, a = a, b = b, S = S_cc, R = R,
                      Ncens = nrow(cens), cens_row = cens[,1], cens_col = cens[,2])

# PRIMARY: censoring above LoD + unstructured residual correlation (data augmentation), log scale
stanmodel_censcorr_A <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens-corr.stan", model_name = "stanmodel_censcorr_A")
samples_censcorr_A <- sampling(stanmodel_censcorr_A, data = standata_cc_A, iter = 20000, warmup = 10000,
                               chains = 5, thin = 10, check_data = FALSE, cores = 5, init_r = 0.5,
                               control = list(adapt_delta = 0.95, max_treedepth = 12))
mcmc_censcorr_A <- do.call(cbind, rstan::extract(samples_censcorr_A,
                                                 pars = c(paste0("gamma[", 1:K, "]"),
                                                          paste0("delta[", 1:K, "]")),
                                                 permuted = TRUE))
save(mcmc_censcorr_A, file = "mcmc/mcmc_censcorr_A.RData")
check_hmc_diagnostics(samples_censcorr_A)
summary(mcmc_censcorr_A)

# SENSITIVITY: log scale, censoring above LoD, DIAGONAL residual covariance
stanmodel_cens_A <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_cens_A")
samples_cens_A <- sampling(stanmodel_cens_A, data = standata_cc_A, iter = 20000, warmup = 10000,
                           chains = 5, thin = 10, check_data = FALSE, cores = 5,
                           control = list(adapt_delta = 0.95, max_treedepth = 12))
mcmc_cens_A <- do.call(cbind, rstan::extract(samples_cens_A,
                                             pars = c(paste0("gamma[", 1:K, "]"),
                                                      paste0("delta[", 1:K, "]")),
                                             permuted = TRUE))
save(mcmc_cens_A, file = "mcmc/mcmc_cens_A.RData")
summary(mcmc_cens_A)

# SENSITIVITY: truncation at LoD (superseded by the censored + correlated primary)
# stanmodel_trunc_A <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_trunc_A")
# samples_trunc_A <- sampling(stanmodel_trunc_A, data = standata_A, iter = 20000, warmup = 10000,
#                             chains = 5, thin = 10, check_data = FALSE, cores = 5,
#                             control = list(adapt_delta = 0.95, max_treedepth = 12))
# mcmc_trunc_A <- do.call(cbind, rstan::extract(samples_trunc_A,
#                                               pars = c(paste0("gamma[", 1:K, "]"),
#                                                        paste0("delta[", 1:K, "]")),
#                                               permuted = TRUE))
# save(mcmc_trunc_A, file = "mcmc/mcmc_trunc_A.RData")
# summary(mcmc_trunc_A)

### Sample B

source("~/Github/ra-biomarker/clean-data-B.R")

# Hyperparameters
a <- b <- rep(0, times = K)
R    <- diag(1e8, nrow = K, ncol = K)   # theta (intercept) prior covariance
S    <- diag(1e8, nrow = K, ncol = K)   # beta/gamma prior covariance (sensitivity)
S_cc <- diag(4, nrow = K, ncol = K)   # weakly-informative slope prior (log scale) for the primary

# Data augmentation indices for the censored cells (D == 1)
cens <- which(D == 1, arr.ind = TRUE)
cens <- cens[order(cens[,1], cens[,2]),]

standata_cc_B <- list(N = N, M = M, K = K, Y = log(Y), U = log(U),
                      t = time, g = diagnosis, id = study_id, a = a, b = b, S = S_cc, R = R,
                      Ncens = nrow(cens), cens_row = cens[,1], cens_col = cens[,2])

# PRIMARY: censoring above LoD + unstructured residual correlation (data augmentation), log scale
stanmodel_censcorr_B <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens-corr.stan", model_name = "stanmodel_censcorr_B")
samples_censcorr_B <- sampling(stanmodel_censcorr_B, data = standata_cc_B, iter = 20000, warmup = 10000,
                               chains = 5, thin = 10, check_data = FALSE, cores = 5, init_r = 0.5,
                               control = list(adapt_delta = 0.95, max_treedepth = 12))
mcmc_censcorr_B <- do.call(cbind, rstan::extract(samples_censcorr_B,
                                                 pars = c(paste0("gamma[", 1:K, "]"),
                                                          paste0("delta[", 1:K, "]")),
                                                 permuted = TRUE))
save(mcmc_censcorr_B, file = "mcmc/mcmc_censcorr_B.RData")
check_hmc_diagnostics(samples_censcorr_B)
summary(mcmc_censcorr_B)

# SENSITIVITY: log scale, censoring above LoD, DIAGONAL residual covariance
stanmodel_cens_B <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-cens.stan", model_name = "stanmodel_cens_B")
samples_cens_B <- sampling(stanmodel_cens_B, data = standata_cc_B, iter = 20000, warmup = 10000,
                           chains = 5, thin = 10, check_data = FALSE, cores = 5,
                           control = list(adapt_delta = 0.95, max_treedepth = 12))
mcmc_cens_B <- do.call(cbind, rstan::extract(samples_cens_B,
                                             pars = c(paste0("gamma[", 1:K, "]"),
                                                      paste0("delta[", 1:K, "]")),
                                             permuted = TRUE))
save(mcmc_cens_B, file = "mcmc/mcmc_cens_B.RData")
summary(mcmc_cens_B)

# SENSITIVITY: truncation at LoD (superseded by the censored + correlated primary)
# stanmodel_trunc_B <- stan_model(file = "~/Github/ra-biomarker/stan/change-point-trunc.stan", model_name = "stanmodel_trunc_B")
# samples_trunc_B <- sampling(stanmodel_trunc_B, data = standata_B, iter = 20000, warmup = 10000,
#                             chains = 5, thin = 10, check_data = FALSE, cores = 5,
#                             control = list(adapt_delta = 0.95, max_treedepth = 12))
# mcmc_trunc_B <- do.call(cbind, rstan::extract(samples_trunc_B,
#                                               pars = c(paste0("gamma[", 1:K, "]"),
#                                                        paste0("delta[", 1:K, "]")),
#                                               permuted = TRUE))
# save(mcmc_trunc_B, file = "mcmc/mcmc_trunc_B.RData")
# summary(mcmc_trunc_B)
