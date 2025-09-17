#-----------------------#
# dataset 2 cleaning
#-----------------------#
library(readxl)
library(tidyverse)
library(Microsoft365R)

file_path <- "~/Documents/RA-Biomarker/data/KevinDat2.xlsx"
data_B_raw <- read_xlsx(file_path)
colnames(data_B_raw) <- tolower(colnames(data_B_raw))

biomarkers <- c("aptivaccp3igg_≥5#00flu", "aptiva_acpafsiggvimentin2_≥5#00au",
                "aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
                "aptivaccp3iga_≥5#00flu", "aptiva_acpafsigavimentin2_≥5#00au",
                "aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
bin_biomarkers <- c("aptivaccp3igg_qual", "aptiva_acpafsiggvimentin2_qual",
                    "aptiva_acpafsiggfibrinogen_qual", "aptiva_acpafsigghistone1_qual", 
                    "aptivaccp3iga_qual", "aptiva_acpafsigavimentin2_qual",
                    "aptiva_acpafsigafibrinogen_qual", "aptiva_acpafsigahistone1_qual")

data_B <- data_B_raw %>% 
  dplyr::rename(diagnosis = casecontrol, 
                study_id = masterstudyid,
                subj_id = masterstudyidnumeric,
                sampnum = sampordernum,
                age = ageserumsample,
                t_days = d_serum_ref) %>% 
  select(study_id, subj_id, sampnum,diagnosis, age, t_days, 
         year_dref, gender, eversmoke, familyhxra, race_ethnic, 
         all_of(biomarkers), all_of(bin_biomarkers))

data_B$subj_id <- ifelse(data_B$diagnosis == "Control", paste0(data_B$subj_id, "_1"), paste(data_B$subj_id)) # control and case had the same subj ID
data_B$subj_id <- as.character(as.integer(factor(data_B$subj_id)))
data_B$study_id <- as.character(as.integer(factor(data_B$study_id)))
data_B$t_yrs <- data_B$t_days/365.25 # changing from days to years

data_B[,biomarkers] <- apply(data_B[,biomarkers], 2, as.numeric)
data_B <- arrange(data_B, as.numeric(subj_id), sampnum) %>% drop_na(all_of(biomarkers))

# age of diagnosis
age_diag_tmp <- data_B %>% 
  arrange(subj_id, t_yrs) %>%
  group_by(subj_id) %>% 
  dplyr::summarize(age_diag = mean(age - t_yrs)) %>%
  ungroup()

data_B <- merge(data_B, age_diag_tmp, by = "subj_id")

rm(age_diag_tmp)

# Outcome matrix
Y <- as.matrix(subset(data_B, select = biomarkers))
colnames(Y) <- c("aptivaccp3igg","aptiva_acpafsiggvimentin2",
                 "aptiva_acpafsiggfibrinogen","aptiva_acpafsigghistone1",
                 "aptivaccp3iga","aptiva_acpafsigavimentin2",
                 "aptiva_acpafsigafibrinogen","aptiva_acpafsigahistone1")
# Y <- apply(Y, 2, function(z) ifelse(z == 0, 1e-6, z))
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(data_B$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
time <- data_B$t_yrs # time before diagnosis
fem <- ifelse(data_B$gender == "F", 1, 0) # indeicator for female
white <- ifelse(data_B$race_ethnic == "W", 1, 0) # indicator for white
black <- ifelse(data_B$race_ethnic == "B", 1, 0) # indicator for black
hispanic <- ifelse(data_B$race_ethnic == "H", 1, 0) # indicator for hispanic
famhx <- ifelse(data_B$familyhxra == "Yes", 1, 0)
age_diag <- data_B$age_diag
diagnosis <- ifelse(data_B$diagnosis == "Case", 1, 0)

# ID vector
subj_id <- as.integer(factor(data_B$subj_id, levels = unique(data_B$subj_id)))
study_id <- as.integer(factor(data_B$study_id, levels = unique(data_B$study_id)))

maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)
D <- apply(Y, 2, function(z) as.numeric(z == max(z)))

# Dichotomus biomarker values

Y_bin <- data_B %>%
  select(all_of(bin_biomarkers))

Y_bin <- as.matrix(Y_bin)
Y_bin <- apply(Y_bin, 2, as.integer)
colnames(Y_bin) <- colnames(Y)
