source ("clean-data.R")

# Add Dichotomous variables

clean$dic_igarfconc_<- ifelse(clean$igarfconc_>8.59151, "pos", "neg")
clean$dic_igmrfconc_<- ifelse(clean$igmrfconc_>26.59495, "pos", "neg")
clean$dic_iggrfconc_<- ifelse(clean$iggrfconc_ >15.79870, "pos", "neg")
clean$dic_igaccpavgconc<-ifelse (clean$igaccpavgconc >110.07650, "pos", "neg")
clean$dic_igmccpavgconc<- ifelse (clean$igmccpavgconc >202.07565, "pos", "neg")
clean$dic_iggccpavgconc<- ifelse (clean$iggccpavgconc >7.66189, "pos", "neg")


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

