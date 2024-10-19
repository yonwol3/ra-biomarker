################################################
##### Simulation study for one outcome ##########
###############################################
library(tidyverse)
sigma_b= 1 # standard deviation for the random intercepts
sigma_e= 2 # standard deviation for the error terms
n=200
m=20
c=10 # num of observations within each subject
kappa=-8 # **** parameter of interest (change point)******
a= -20 # parameter for the uniform distribution 
b= 5 # parameter for the uniform distribution
b_1= 0.1 # time parameter
b_2= 0.6 # trt paramter 
b_3= 0.3 # ****parameter of interest (interaction term)*****
# Ei ~ N(0,sigma_e^2)
# alphai ~ N(0,sigma_b^2)
# timeij ~ unif (a,b)
# We have 20 subject IDs and 1:1 ratio for treatment allocation
observation<-1:200
id<-1:20 
time<-matrix(NA, nrow = m, ncol=c)
for (i in 1:length(id)) {
  time[i,]<-runif(n=c,a,b)
}
time<-as.vector(t(time))
df<- data.frame (observation=observation, id=rep(id,each=c),
                 time=time)
df<- df %>% mutate (ccp3=ifelse(id %in% c(1:10),"1","0"))
df$alpha=rep(rnorm(m, mean=0, sd=sigma_b), each=c)
df$ccp3<-as.numeric(df$ccp3)
iga=vector() # outcome 
for (i in 1:n) {
  if (df[i,3]<=kappa | df[i,4]==0) {
        
            iga[i]<-df[i,5] + b_1*df[i,3] + b_2*df[i,4] + rnorm(1,mean = 0, sd=sigma_e)
      
            } else {
    
          iga[i]<-df[i,5] + b_1*df[i,3] + b_2*df[i,4] + b_3*(df[i,4]-kappa) +rnorm(1,mean = 0, sd=sigma_e)
    }
}

df$outcome=iga


### Now create a simulation loop and save the mean, median, zero contained or not,
# true_value is contained or not for each simulation for the two parameters of interest
# so use two different data_frames to save the data 
# then calculate bias, mse, power, coverage probability for each parameter (b_3 and kappa)

