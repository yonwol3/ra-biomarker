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
data_A <- read.delim("~/Documents/RA-Biomarker/data/forKevin.txt", stringsAsFactors = FALSE)
names(data_A) <- tolower(names(data_A))

serum_names <- c("igarfconc_", "igmrfconc_", "iggrfconc_",
                 "igaccpavgconc", "igmccpavgconc", "iggccpavgconc")

# Fix Race Variable

data_A_case <- subset(data_A, diagnosis == "RA")
data_A_case$race_ethnic <- NULL
data_A_control <- subset(data_A, diagnosis == "Control")

data_A_case <- data_A_control %>% 
  group_by(study_id) %>%
  filter(row_number() == 1) %>%
  dplyr::summarize(race_ethnic = race_ethnic) %>%
  ungroup() %>%
  inner_join(data_A_case, by = "study_id")

data_A <- rbind(data_A_case, data_A_control)

rm(data_A_case, data_A_control)

## Generate Data Frame

data_A <- data_A[order(data_A$subj_id, data_A$sampnum),]
data_A <- data_A[complete.cases(data_A[,serum_names]),] # remove one observation missing 3 measurements
data_A$subj_id <- as.integer(factor(data_A$subj_id))
data_A$study_id <- as.integer(factor(data_A$study_id))
data_A$diagnosis[is.na(data_A$diagnosis)] <- "RA"

# Getting age at diagnosis variable from age and t_years variables
age_diag_tmp <- data.frame(age_diag = data_A$age - data_A$t_yrs, subj_id = data_A$subj_id)
mean_age_diag <- aggregate(age_diag ~ subj_id, age_diag_tmp, mean)
data_A <- join(data_A, mean_age_diag, by = "subj_id", type = "left", match = "all")

rm(age_diag_tmp, mean_age_diag)

## Generate objects for data list to be passed to JAGS

# Outcome matrix
Y <- as.matrix(subset(data_A, select = serum_names) )
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(data_A$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
time <- data_A$t_yrs # time before diagnosis
fem <- ifelse(data_A$gender == "F", 1, 0) # indicator for female
white <- ifelse(data_A$race_ethnic == "W", 1, 0) # indicator for white
black <- ifelse(data_A$race_ethnic == "B", 1, 0) # indicator for black
hisp <- ifelse(data_A$race_ethnic == "H", 1, 0) # indicator for Hispanic
famhx <- ifelse(data_A$familyhxra == "Yes", 1, 0)
age_diag <- data_A$age_diag
diagnosis <- ifelse(data_A$diagnosis == "RA", 1, 0)

# ID vector
subj_id <- data_A$subj_id
study_id <- data_A$study_id

maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)
D <- apply(Y, 2, function(z) as.numeric(z == max(z)))

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
