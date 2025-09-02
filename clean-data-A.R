#############################################
## PURPOSD: Script for cleaning data       ##
##          in RA serum sample analysis    ##
## BY:      Kevin Josey                    ##
#############################################

require(tidyverse)
require(plyr)
library(readxl)

## Data Cleaning

setwd("~/Documents/RA-Biomarker/")
raDat <- read.delim("~/Documents/RA-Biomarker/data/forKevin.txt", stringsAsFactors = FALSE)
names(raDat) <- tolower(names(raDat))

# raDat$race_ethnic<-race_clean$RACE_ETHNIC # the correct race/ethnicity values
raDat_case <- subset(raDat, diagnosis == "RA")
raDat_control <- subset(raDat, diagnosis == "Control")

serum_names <- c("igarfconc_", "igmrfconc_", "iggrfconc_",
                 "igaccpavgconc", "igmccpavgconc", "iggccpavgconc")

raDat <- raDat[order(raDat$subj_id, raDat$sampnum),]
raDat <- raDat[complete.cases(raDat[,serum_names]),] # remove one observation missing 3 measurements

raDat$subj_id <- as.integer(factor(raDat$subj_id))
raDat$study_id <- as.integer(factor(raDat$study_id))

## Generate Data Frame for Table 1

# Getting age at diagnosis variable from age and t_years variables
tmp <- data.frame(agediag = raDat$age - raDat$t_yrs, subj_id = raDat$subj_id)
mean_agediag <- aggregate(agediag ~ subj_id, tmp, mean)
raDat <- join(raDat, mean_agediag, by = "subj_id", type = "left", match = "all")

## Generate objects for data list to be passed to JAGS

# Outcome matrix
Y <- as.matrix(subset(raDat, select = serum_names) )
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(raDat$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
time <- raDat$t_yrs # time before diagnosis
fem <- ifelse(raDat$gender == "F", 1, 0) # indicator for female
white <- ifelse(raDat$race_ethnic == "W", 1, 0) # indicator for white
black <- ifelse(raDat$race_ethnic == "B", 1, 0) # indicator for black
hisp <- ifelse(raDat$race_ethnic == "H", 1, 0) # indicator for Hispanic
famhx <- ifelse(raDat$familyhxra == "Yes", 1, 0)

age_diag <- raDat$agediag
diagnosis <- ifelse(raDat$diagnosis == "RA", 1, 0)

rm(tmp, mean_agediag)

# ID vector
subj_id <- raDat$subj_id
study_id <- raDat$study_id

maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)
D <- apply(Y, 2, function(z) as.numeric(z == max(z)))

data_A <- data.frame(time, age_diag, fem, white, black, hisp, famhx, subj_id, study_id, diagnosis, Y)

## Add Dichotomous variables

data_A$bin_igarfconc_ <- ifelse(data_A$igarfconc_ > 8.59151, 1, 0)
data_A$bin_igmrfconc_ <- ifelse(data_A$igmrfconc_ > 26.59495, 1, 0)
data_A$bin_iggrfconc_ <- ifelse(data_A$iggrfconc_ > 15.79870, 1, 0)
data_A$bin_igaccpavgconc <- ifelse (data_A$igaccpavgconc > 110.07650, 1, 0)
data_A$bin_igmccpavgconc <- ifelse (data_A$igmccpavgconc > 202.07565, 1, 0)
data_A$bin_iggccpavgconc <- ifelse (data_A$iggccpavgconc > 7.66189, 1, 0)

# data_A_sens <- data_A %>% filter(time<=0)

Y_bin <- data_A %>% 
  select(bin_igarfconc_, bin_igmrfconc_, bin_iggrfconc_, 
         bin_igaccpavgconc, bin_igmccpavgconc, bin_iggccpavgconc) 

Y_bin <- as.matrix(Y_bin)
Y_bin <- apply(Y_bin, 2, as.integer)
colnames(Y_bin) <- serum_names 
