#############################################
## PURPOSD: Script for cleaning data       ##
##          in RA serum sample analysis    ##
## BY:      Kevin Josey                    ##
#############################################

require(plyr)
require(tidyverse)

## Data Cleaning
setwd("~/Dropbox (Personal)/Projects/RA-Biomarker/")

raDat_raw <- read.delim("data/forKevin.txt", stringsAsFactors = FALSE)
names(raDat_raw) <- tolower(names(raDat_raw))

raDat_case <- subset(raDat_raw, diagnosis == "RA")
raDat_control <- subset(raDat_raw, diagnosis == "Control")

serum_names <- names(raDat_case)[8:13]

# Wide to Long Format
raDat_race <- subset(raDat_control, select = c(study_id, race_ethnic))
raDat_tmp1 <- raDat_race[!duplicated(raDat_race$study_id),]
raDat_tmp2 <- subset(raDat_case, select = -c(race_ethnic, subj_id))

# 5 Patients Missing Race
raDat <- merge(raDat_tmp1, raDat_tmp2, by = "study_id", all.y = TRUE) 

rm(raDat_raw, raDat_race, raDat_tmp1, raDat_tmp2) # remove temporary data frames

# raDat <- raDat[!is.na(raDat$race),] # remove responses missing race
# raDat <- raDat[!is.na(raDat$familyhxra),] # remove responses missing family history
# apply(raDat, 2, function(x) sum(is.na(x)) )

raDat <- raDat[order(raDat$study_id, raDat$sampnum),]
raDat <- raDat[complete.cases(raDat[,serum_names]),] # remove one observation missing 3 measurements

raDat$study_id <- as.integer(factor(raDat$study_id))

## Generate Data Frame for Table 1

# Getting age at diagnosis variable from age and t_years variables
raDat$age_diag <- raDat$age - raDat$t_yrs # slight differences due to rounding
mean_agediag <- aggregate(age_diag ~ study_id, raDat, mean)
colnames(mean_agediag) <- c("study_id", "mean_agediag")
raDat_case_join <- join(raDat, mean_agediag, by = "study_id", type = "left", match = "all")

# Centering age at diagnosis variable
raDat$agediag_c <- raDat$age_diag - mean(raDat$age_diag)

# Creating a white/not white race/ethnicity variable
raDat_case_join$white <- ifelse(raDat_case_join$race_ethnic == "W", 1, 0)

# Deleting all but first observation per subject
raDat_t1 <- raDat_case_join[match(unique(raDat_case_join$study_id), raDat_case_join$study_id),]

rm(mean_agediag, raDat_case_join)

## Generate objects for data list to be passed to JAGS

# Outcome matrix
Y <- as.matrix( subset(raDat, select = c(igarfconc_, igmrfconc_, iggrfconc_,
                                         igaccpavgconc, igmccpavgconc, iggccpavgconc)) )
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(raDat$study_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
t <- raDat$t_yrs # time before diagnosis
fem <- ifelse(raDat$gender == "F", 1, 0) # indeicator for female
nw <- ifelse(raDat$race_ethnic == "W", 0, 1) # indicator for non-white
famhx <- ifelse(raDat$familyhxra == "No", 0, 1)
bage.tmp <- raDat$agediag_c[unique(raDat$study_id)]
bage <- rep(bage.tmp, times = table(raDat$study_id))

rm(bage.tmp)

# ID vector
study_id <- raDat$study_id
clean <- data.frame(t, bage, study_id, Y)
write.csv(clean, "data/clean.csv")
