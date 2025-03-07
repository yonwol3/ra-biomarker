#############################################
## PURPOSD: Script for cleaning data       ##
##          in RA serum sample analysis    ##
## BY:      Kevin Josey                    ##
#############################################

require(plyr)
require(tidyverse)

## Data Cleaning
setwd("~/Documents/RA-Biomarker/")

raDat <- read.delim("data/forKevin.txt", stringsAsFactors = FALSE)
names(raDat) <- tolower(names(raDat))

raDat_case <- subset(raDat, diagnosis == "RA")
raDat_control <- subset(raDat, diagnosis == "Control")

serum_names <- names(raDat_case)[8:13]

raDat <- raDat[order(raDat$subj_id, raDat$sampnum),]
raDat <- raDat[complete.cases(raDat[,serum_names]),] # remove one observation missing 3 measurements

raDat$subj_id <- as.integer(factor(raDat$subj_id))

## Generate Data Frame for Table 1

# Getting age at diagnosis variable from age and t_years variables
tmp <- data.frame(agediag = raDat$age - raDat$t_yrs, subj_id = raDat$subj_id)
mean_agediag <- aggregate(agediag ~ subj_id, tmp, mean)
raDat <- join(raDat, mean_agediag, by = "subj_id", type = "left", match = "all")

## Generate objects for data list to be passed to JAGS

# Outcome matrix
Y <- as.matrix( subset(raDat, select = serum_names) )
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(raDat$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
time <- raDat$t_yrs # time before diagnosis
fem <- ifelse(raDat$gender == "F", 1, 0) # indeicator for female
white <- ifelse(raDat$race_ethnic == "W", 1, 0) # indicator for white
black <- ifelse(raDat$race_ethnic == "B", 1, 0) # indicator for black
hisp <- ifelse(raDat$race_ethnic == "H", 1, 0) # indicator for hispanic
famhx <- ifelse(raDat$familyhxra == "Yes", 1, 0)
age_diag <- raDat$agediag
diagnosis <- ifelse(raDat$diagnosis == "RA", 1, 0)

rm(tmp, mean_agediag, bage.tmp)

# ID vector
subj_id <- raDat$subj_id
study_id <- raDat$study_id

cens_max <- apply(Y, 2, function(z) as.numeric(z == max(z)))
cens_min <- apply(Y, 2, function(z) as.numeric(z == min(z)))
maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)

