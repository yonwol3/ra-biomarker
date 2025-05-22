#############################################
## PURPOSD: Script for cleaning data       ##
##          in RA serum sample analysis    ##
## BY:      Kevin Josey                    ##
#############################################

library(Microsoft365R)
require(tidyverse)
require(plyr)
library(readxl)

## Data Cleaning

# onedrive <- get_business_onedrive()
# file_path <- "Attachments/forKevin.txt "
# file_path("~/Documents/RA-Biomarker/")
# temp_file <- tempfile(fileext = ".txt")
# onedrive$download_file(
#   src = file_path,
#   dest = temp_file,
#   overwrite = TRUE
# )
# raDat <- read.delim(temp_file, stringsAsFactors = FALSE)
# names(raDat) <- tolower(names(raDat))
# unlink(temp_file)

# onedrive <- get_business_onedrive()
# file_path <- "Attachments/DOD_AddRace.xls"
# temp_file <- tempfile(fileext = ".xls")
# onedrive$download_file(
#   src = file_path,
#   dest = temp_file,
#   overwrite = TRUE
# )
# race_clean <- read_xls(temp_file, sheet=2)
# unlink(temp_file)

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

clean <- data.frame(time, age_diag, fem, white, black, hisp, famhx, subj_id, study_id, diagnosis, Y)

## Add Dichotomous variables

clean$bin_igarfconc_ <- ifelse(clean$igarfconc_ > 8.59151, 1, 0)
clean$bin_igmrfconc_ <- ifelse(clean$igmrfconc_ > 26.59495, 1, 0)
clean$bin_iggrfconc_ <- ifelse(clean$iggrfconc_ > 15.79870, 1, 0)
clean$bin_igaccpavgconc <- ifelse (clean$igaccpavgconc > 110.07650, 1, 0)
clean$bin_igmccpavgconc <- ifelse (clean$igmccpavgconc > 202.07565, 1, 0)
clean$bin_iggccpavgconc <- ifelse (clean$iggccpavgconc > 7.66189, 1, 0)

# clean_sens<- clean %>% filter(time<=0)

Y_bin <- clean %>% 
  select(bin_igarfconc_, bin_igmrfconc_, bin_iggrfconc_, 
         bin_igaccpavgconc, bin_igmccpavgconc, bin_iggccpavgconc) 

Y_bin <- as.matrix(Y_bin)
Y_bin <- apply(Y_bin, 2, as.integer)
colnames(Y_bin) <- serum_names 
