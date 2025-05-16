source ("~/Github/ra-biomarker/clean-data.R")

# Add Dichotomous variables

clean$dic_igarfconc_ <- ifelse(clean$igarfconc_ > 8.59151, 1, 0)
clean$dic_igmrfconc_ <- ifelse(clean$igmrfconc_ > 26.59495, 1, 0)
clean$dic_iggrfconc_ <- ifelse(clean$iggrfconc_ > 15.79870, 1, 0)
clean$dic_igaccpavgconc <- ifelse (clean$igaccpavgconc > 110.07650, 1, 0)
clean$dic_igmccpavgconc <- ifelse (clean$igmccpavgconc > 202.07565, 1, 0)
clean$dic_iggccpavgconc <- ifelse (clean$iggccpavgconc > 7.66189, 1, 0)

# remove observations that are t>0 (after diagnosis)
# remove from both cases and controls

# clean_sens<- clean %>% filter(time<=0)

Y <- clean %>% 
  select(dic_igarfconc_, dic_igmrfconc_, dic_iggrfconc_, 
         dic_igaccpavgconc, dic_igmccpavgconc, dic_iggccpavgconc) %>% 
  as.matrix()

cens_max <- apply(Y, 2, function(z) as.numeric(z == max(z)))
cens_min <- apply(Y, 2, function(z) as.numeric(z == min(z)))
maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)

