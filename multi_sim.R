##########################################
## multivariate simulation with stan
## change_point a (correlation on error)
## By: Yonatan Wolde
##########################################

library(tidyverse)
library(MASS)
library(rstan)

set.seed(2222)

#-------------------------#
# Parameter Specification
#-------------------------#


a <- -20      # Uniform distribution lower bound
b <- 5        # Uniform distribution upper bound
n <- 200      # Total number of observations
m <- 50       # Number of subjects
c <- n/m     # Observations per subject
k<-6         # of biomarkers
sigma_0 <- 1  # Standard deviation for random intercepts
sigma_e <- 2  # Standard deviation for error terms
var_cov_sigma_0<-diag(sigma_0^2,k,k) # var-cov for alpha
var_cov_sigma_e<- matrix(rexp(k*k, rate = 4), nrow = k)
diag(var_cov_sigma_e)<-rep(sigma_e^2,k) # correlation of the errors
beta1<- rnorm(k)  # Time parameters
beta2<- rnorm(k)  # Treatment parameters
gamma<- rnorm(k) # Interaction parameters
kappa <- runif(6, min=a, max=b)  # Change point parameters
mu<- rnorm(k)


#---------------------------#
# Simulation study 
#--------------------------# 

Y<- matrix (NA, nrow=n, ncol = k)
nsim= 1  

  id <- rep(1:m, each=c)
  time <- runif(n, a, b)
  ccp3 <- as.numeric(id <= 20)
  alpha<- mvrnorm(m, mu = mu, Sigma = var_cov_sigma_0) # random intercept 
  error<-mvrnorm(n, mu=rep(0,k), Sigma= var_cov_sigma_e) # random errors
  
  for (j in 1:k) {
    Y[,j]<- rep(alpha[, j],each=c) + beta1[j]*time+beta2[j]*ccp3 + 
             ifelse(time > kappa[j] & ccp3 == 1, gamma[j] * (time - kappa[j]), 0) + error[j]
    
  }
  
  dat<- as.data.frame(cbind(Y,id=id, time=time, ccp3=ccp3))
  
  stan_dat<- list(N=nrow(dat),
                  M=length(unique(dat$id)),
                  K=k, 
                  id=dat$id,
                  t=dat$time,
                  g=dat$ccp3,
                  a=0, 
                  b=0,
                  R= diag(10^6, nrow= k),
                  S=diag(10^6, nrow = k),
                  Y= Y)
  
                  
  stan_fit<- stan(data = stan_dat, file = "Stan/multi_sim.stan", 
                  chains = 1, iter = 4000, warmup = 2000)
  
  
  

