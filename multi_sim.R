##########################################
## multivariate simulation with stan
## change_point a (correlation on error)
## By: Yonatan Wolde
##########################################

library(tidyverse)
library(MASS)

set.seed(2222)

#-------------------------#
# Parameter Specification
#-------------------------#


a <- -20      # Uniform distribution lower bound
b <- 5        # Uniform distribution upper bound
sigma_0 <- 1  # Standard deviation for random intercepts
sigma_e <- 2  # Standard deviation for error terms
var_cov_sigma_0<-diag(sigma_0^2,k,k) # var-cov for alpha
var_cov_sigma_e<- matrix(rexp(k*k, rate = 4), nrow = k)
diag(var_cov_sigma_e)<-rep(sigma_e^2,k) # correlation of the errors
n <- 200      # Total number of observations
m <- 50       # Number of subjects
c <- n/m     # Observations per subject
k<-6      # of biomarkers
beta1<- rnorm(k)  # Time parameters
beta2<- rnorm(k)  # Treatment parameters
gamma<- rnorm(k) # Interaction parameters
kappa <- runif(6, min=a, max=b)  # Change point parameter
mu<- rnorm(k)


#---------------------------#
# Simulation study
#--------------------------# 
nsim= 1  
for (i in 1:nsim) {
  id <- rep(1:m, each=c)
  time <- runif(n, a, b)
  ccp3 <- as.numeric(id <= 20)
  alpha<- mvrnorm(m, mu = mu, Sigma = var_cov_sigma_0) # random intercept 
  error<-mvrnorm(n, mu=rep(0,k), Sigma= var_cov_sigma_e) # random errors
  
  for (j in 1:k) {
    y[j]<- rep(alpha[, j],each=c) + beta1[j]*time+beta2[j]*ccp3 + 
             ifelse(time > kappa[j] & ccp3 == 1, gamma[j] * (time - kappa[j]), 0) + error[j]
      
  }
  
  y_1<- alpha
  
}
