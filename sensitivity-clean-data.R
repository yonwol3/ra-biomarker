source ("clean-data.R")

# remove observations that are t>0 (after diagnosis)
# remove from both cases and controls

clean_sens<- clean %>% 
              filter(time<=0)

Y <- as.matrix(clean_sens[, 11:16])
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(clean_sens$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

cens_max <- apply(Y, 2, function(z) as.numeric(z == max(z)))
cens_min <- apply(Y, 2, function(z) as.numeric(z == min(z)))
maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)

