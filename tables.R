#############################################
## PURPOSE: Script for creation of tables  ##
##          in RA data analysis            ##
## BY:      Yonatan Wolde                 ##
#############################################

# Libraries

library(plyr)
library(tidyverse)
library(labelled)
library(table1)
library(naniar)
library(mcmcse)

# Set working directory
#setwd("~/Dropbox/Projects/RA-Biomarker/")

hpd = function(x, alpha = 0.05) {
  n = length(x)
  m = round(n * alpha)
  x = sort(x)
  y = x[(n - m + 1):n] - x[1:m]
  z = min(y)
  k = which(y == z)[1]
  c(x[k], x[n - m + k])
}

# Load the dataframes
source("clean-data.R")
source("clean-data-new.R")

clean$agediag<-raDat$agediag # adding mean age at diagnosis from raDat
dat_2$agediag<-dat_2$age-dat_2$t_yrs # age at diagnosis new dataset

# covariate elements need some renaming

clean$fem=factor(clean$fem,labels =c("Male","Female"))
clean$nw=factor(clean$nw,labels =c("White","Non-white"))
clean$eversmoke=raDat$eversmoke # smoking status from raDat
clean$diagnosis= factor(clean$diagnosis,labels=c("No RA", "RA"))
clean <- clean %>%
  dplyr::mutate(
    famhx_c = case_when(
      famhx == 0 ~ "No",
      famhx == 1 ~ "Yes",
      is.na(famhx) ~ "Not available"
    ), 
    smoking= case_when(
      eversmoke == "No"~ "No",
      eversmoke == "Yes"~ "Yes",
      TRUE ~ "Not available"  
    )
  )

dat_2$gender=factor(dat_2$gender,levels=c("M","F"),labels =c("Male","Female"))
dat_2$nw=ifelse(dat_2$race_ethnic=="W",0,1)# non-white indicator
dat_2$nw=factor(dat_2$nw,labels =c("White","Non-white")) 
dat_2$diagnosis=factor(dat_2$diagnosis, levels=c("Control","Case"), labels=c("No RA", "RA"))

dat_2 <- dat_2 %>%
  dplyr::mutate(
    famhx_c = case_when(
      familyhxra == "No" ~ "No",
      familyhxra == "Yes"~ "Yes",
      TRUE ~ "Not available"
    ), 
    smoking= case_when(
      eversmoke == "No"~ "No",
      eversmoke == "Yes"~ "Yes",
      TRUE ~ "Not available"  
    )
  )

# Take mean values of biomarkes per subject to create table-1

clean_tb1 <- clean %>% 
  group_by(subj_id) %>%
  dplyr::summarise(
    gender = first(fem),
    race = first(nw),
    famhx_c = first(famhx_c),
    diagnosis=first(diagnosis),
    smoking=first(smoking),
    mean_agediag = mean(agediag),
    across(starts_with("ig"), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}"),
  )

dat_2_tb1<-dat_2 %>% 
  group_by(subj_id) %>%
  dplyr::summarise(
    gender = first(gender),
    race = first(nw),
    famhx_c = first(famhx_c),
    smoking=first(smoking),
    diagnosis=first(diagnosis),
    mean_agediag = mean(agediag),
    across(starts_with("apti"), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}"),
  )


# specifying variable names shown in the final table-1 output

clean_tb1 <- set_variable_labels(clean_tb1,
                           gender= "Sex",
                           race= "Race",
                           famhx_c = "Family RA History",
                           smoking= "Eversmoke", 
                           mean_agediag="Age at Diagnosis",
                           mean_igarfconc_="RF IgA",
                           mean_igmrfconc_="RF IgM",
                           mean_iggrfconc_="RF IgG",
                           mean_igaccpavgconc="ACPA IgA",
                           mean_igmccpavgconc="ACPA IgM",
                           mean_iggccpavgconc="ACPA IgG")


# Then create your Table 1
covariates <- c("gender", "race","famhx_c","smoking")
biomarkers <- c("mean_agediag","mean_igarfconc_", "mean_igmrfconc_", "mean_iggrfconc_", 
                "mean_igaccpavgconc","mean_igmccpavgconc","mean_iggccpavgconc" )
t1_variables<-c(covariates,biomarkers)
categorical_vars <-covariates


# Enclose variable names with backticks to handle spaces
t1_variables_backticked <- paste0("`", t1_variables, "`")

# Create the table with the corrected formula
t1 <- table1(
  as.formula(paste0("~", paste0(t1_variables_backticked, collapse = "+"), "|", "diagnosis")),
  data = clean_tb1,
  caption = "Table1: Descriptive summary of data by RA status"
)

t1  # Display the table





## Table 

load("mcmc/mcmc_1b.RData")
load("mcmc/mcmc_2.RData")
load("mcmc/mcmc_3b.RData")
load("mcmc/mcmc_4.RData")

kappa.names <- c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]", "kappa[5]", "kappa[6]")

kappa_1 <- mcmc_1b[[1]][,kappa.names]
kappa_2 <- mcmc_2[[1]][,kappa.names]
kappa_3 <- mcmc_3b[[1]][,kappa.names]
kappa_4 <- mcmc_4[[1]][,kappa.names]

# Posterior Probabilities
idx <- 1:6
grid.tmp <- expand.grid(idx, idx)
grid.idx <- grid.tmp[grid.tmp[,1] < grid.tmp[,2],]
post.prob.mat <- matrix("", nrow = nrow(grid.idx), ncol = 5)
colnames(post.prob.mat) <- c("Comparison", "Model (1)", "Model (2)", "Model (3)", "Model (4)")

for (k in 1:nrow(grid.idx)) {
  
  post.prob.mat[k, 1] <- paste("Pr{", grid.idx[k,1], "<", grid.idx[k,2], "|y}", sep = "")
  post.prob.mat[k, 2] <- mean(as.numeric(kappa_1[,grid.idx[k,1]] < kappa_1[,grid.idx[k,2]]))
  post.prob.mat[k, 3] <- mean(as.numeric(kappa_2[,grid.idx[k,1]] < kappa_2[,grid.idx[k,2]]))
  post.prob.mat[k, 4] <- mean(as.numeric(kappa_3[,grid.idx[k,1]] < kappa_3[,grid.idx[k,2]]))
  post.prob.mat[k, 5] <- mean(as.numeric(kappa_4[,grid.idx[k,1]] < kappa_4[,grid.idx[k,2]]))
  
}

write.csv(post.prob.mat, file = "data/table-2.csv")

## Table 3

kappa.tbl <- matrix("", nrow = 6, ncol = 20)

kmcse.tmp1 <- apply(kappa_1, 2, mcse)
kmcse.tmp2 <- apply(kappa_2, 2, mcse)
kmcse.tmp3 <- apply(kappa_3, 2, mcse)
kmcse.tmp4 <- apply(kappa_4, 2, mcse)

for(j in 1:nrow(kappa.tbl)) {
  
  kappa.tbl[j,2] <- as.character(round(kmcse.tmp1[[j]][[1]], 3))
  kappa.tbl[j,7] <- as.character(round(kmcse.tmp2[[j]][[1]], 3))
  kappa.tbl[j,12] <- as.character(round(kmcse.tmp3[[j]][[1]], 3))
  kappa.tbl[j,17] <- as.character(round(kmcse.tmp4[[j]][[1]], 3))
  
  kappa.tbl[j,3] <- as.character(round(kmcse.tmp1[[j]][[2]], 3))
  kappa.tbl[j,8] <- as.character(round(kmcse.tmp2[[j]][[2]], 3))
  kappa.tbl[j,13] <- as.character(round(kmcse.tmp3[[j]][[2]], 3))
  kappa.tbl[j,18] <- as.character(round(kmcse.tmp4[[j]][[2]], 3))
  
}

kappa.tbl[,4] <- as.character(round(apply(kappa_1, 2, sd), 3))
kappa.tbl[,9] <- as.character(round(apply(kappa_2, 2, sd), 3))
kappa.tbl[,14] <- as.character(round(apply(kappa_3, 2, sd), 3))
kappa.tbl[,19] <- as.character(round(apply(kappa_4, 2, sd), 3))

khpd.int1 <- t(apply(kappa_1, 2, hpd))
khpd.int2 <- t(apply(kappa_2, 2, hpd))
khpd.int3 <- t(apply(kappa_3, 2, hpd))
khpd.int4 <- t(apply(kappa_4, 2, hpd))

for (i in 1:nrow(kappa.tbl)) {
  
  kappa.tbl[i,5] <- paste("(", round(khpd.int1[i,1], 3), ", ", 
                       round(khpd.int1[i,2], 3), ")", sep = "")
  
  kappa.tbl[i,10] <- paste("(", round(khpd.int2[i,1], 3), ", ", 
                        round(khpd.int2[i,2], 3), ")", sep = "")
  
  kappa.tbl[i,15] <- paste("(", round(khpd.int3[i,1], 3), ", ", 
                        round(khpd.int3[i,2], 3), ")", sep = "")
  
  kappa.tbl[i,20] <- paste("(", round(khpd.int4[i,1], 3), ", ", 
                        round(khpd.int4[i,2], 3), ")", sep = "")
  
}

kappa.tbl[,1] <- colnames(kappa_1)
kappa.tbl[,6] <- colnames(kappa_2)
kappa.tbl[,11] <- colnames(kappa_3)
kappa.tbl[,16] <- colnames(kappa_4)

colnames(kappa.tbl) <- rep(c("Parameter", "Posterior Mean", "MC Std. Err.", "Std. Dev.", "95% HPDI"), 4)

write.csv(kappa.tbl, file = "data/table-3.csv")

## Table 4

fix_eff_1 <- mcmc_1b[[1]][,-which(colnames(mcmc_1[[1]]) %in% kappa.names)]
fix_eff_2 <- mcmc_2[[1]][,-which(colnames(mcmc_2[[1]]) %in% kappa.names)]
fix_eff_3 <- mcmc_3b[[1]][,-which(colnames(mcmc_3[[1]]) %in% kappa.names)]
fix_eff_4 <- mcmc_4[[1]][,-which(colnames(mcmc_4[[1]]) %in% kappa.names)]

# Fixed Effects Table

fe.tbl <- matrix("", nrow = ncol(fix_eff_1), ncol = 20)

fix_eff_1 <- fix_eff_1[,c(31, 1:5, 32, 6:10, 33, 11:15, 34, 16:20,
                          35, 21:25, 36, 26:30)]
fix_eff_2 <- fix_eff_2[,c(19:21, 1:3, 22:24, 4:6, 25:27, 7:9, 28:30,
                          10:12, 31:33, 13:15, 34:36, 16:18)]
fix_eff_3 <- fix_eff_3[,c(31, 1:5, 32, 6:10, 33, 11:15, 34, 16:20,
                          35, 21:25, 36, 26:30)]
fix_eff_4 <- fix_eff_4[,c(19:21, 1:3, 22:24, 4:6, 25:27, 7:9, 28:30,
                          10:12, 31:33, 13:15, 34:36, 16:18)]

mcse.tmp1 <- apply(fix_eff_1, 2, mcse)
mcse.tmp2 <- apply(fix_eff_2, 2, mcse)
mcse.tmp3 <- apply(fix_eff_3, 2, mcse)
mcse.tmp4 <- apply(fix_eff_4, 2, mcse)

for(j in 1:nrow(fe.tbl)) {
  
  fe.tbl[j,2] <- as.character(round(mcse.tmp1[[j]][[1]], 3))
  fe.tbl[j,7] <- as.character(round(mcse.tmp2[[j]][[1]], 3))
  fe.tbl[j,12] <- as.character(round(mcse.tmp3[[j]][[1]], 3))
  fe.tbl[j,17] <- as.character(round(mcse.tmp4[[j]][[1]], 3))
  
  fe.tbl[j,3] <- as.character(round(mcse.tmp1[[j]][[2]], 3))
  fe.tbl[j,8] <- as.character(round(mcse.tmp2[[j]][[2]], 3))
  fe.tbl[j,13] <- as.character(round(mcse.tmp3[[j]][[2]], 3))
  fe.tbl[j,18] <- as.character(round(mcse.tmp4[[j]][[2]], 3))
  
}

fe.tbl[,4] <- as.character(round(apply(fix_eff_1, 2, sd), 3))
fe.tbl[,9] <- as.character(round(apply(fix_eff_2, 2, sd), 3))
fe.tbl[,14] <- as.character(round(apply(fix_eff_3, 2, sd), 3))
fe.tbl[,19] <- as.character(round(apply(fix_eff_4, 2, sd), 3))

hpd.int1 <- t(apply(fix_eff_1, 2, hpd))
hpd.int2 <- t(apply(fix_eff_2, 2, hpd))
hpd.int3 <- t(apply(fix_eff_3, 2, hpd))
hpd.int4 <- t(apply(fix_eff_4, 2, hpd))

for (i in 1:nrow(fe.tbl)) {
  
  fe.tbl[i,5] <- paste("(", round(hpd.int1[i,1], 3), ", ", 
                       round(hpd.int1[i,2], 3), ")", sep = "")
  
  fe.tbl[i,10] <- paste("(", round(hpd.int2[i,1], 3), ", ", 
                        round(hpd.int2[i,2], 3), ")", sep = "")
  
  fe.tbl[i,15] <- paste("(", round(hpd.int3[i,1], 3), ", ", 
                        round(hpd.int3[i,2], 3), ")", sep = "")
  
  fe.tbl[i,20] <- paste("(", round(hpd.int4[i,1], 3), ", ", 
                        round(hpd.int4[i,2], 3), ")", sep = "")
  
}

fe.tbl[,1] <- colnames(fix_eff_1)
fe.tbl[,6] <- colnames(fix_eff_2)
fe.tbl[,11] <- colnames(fix_eff_3)
fe.tbl[,16] <- colnames(fix_eff_4)

colnames(fe.tbl) <- rep(c("Parameter", "Posterior Mean", "MC Std. Err.", "Std. Dev.", "95% HPDI"), 4)

write.csv(fe.tbl, file = "data/table-4.csv")
