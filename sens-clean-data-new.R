source ("clean-data-new.R")

# remove observations that are t>0 (after diagnosis)
# remove from both cases and controls

dat2_sens<- dat_2 %>% 
  filter(t_yrs<=0)

Y <- as.matrix(dat2_sens[, 12:23])
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(dat2_sens$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

cens_max <- apply(Y, 2, function(z) as.numeric(z == max(z)))
cens_min <- apply(Y, 2, function(z) as.numeric(z == min(z)))
maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)