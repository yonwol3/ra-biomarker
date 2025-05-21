#-----------------------#
# dataset 2 cleaning
#-----------------------#
library(readxl)
library(tidyverse)
library(Microsoft365R)

# onedrive<- get_business_onedrive()
# file_path <- "Attachments/KevinDat2.xlsx"
# temp_file <- tempfile(fileext = ".xlsx")
# onedrive$download_file(
#   src = file_path,
#   dest = temp_file,
#   overwrite = TRUE
# )
# dat_2raw <- read_xlsx(temp_file)
# unlink(temp_file)

file_path <- "~/Documents/RA-Biomarker/data/KevinDat2.xlsx"
dat_2_raw <- read_xlsx(file_path)
colnames(dat_2_raw) <- tolower(colnames(dat_2_raw))

biomarkers<-c("aptivaccp3igg_≥5#00flu", "aptiva_acpafsiggvimentin2_≥5#00au",
              "aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
              "aptivaccp3iga_≥5#00flu", "aptiva_acpafsigavimentin2_≥5#00au",
              "aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
bin_biomarkers <- c("aptivaccp3igg_qual", "aptiva_acpafsiggvimentin2_qual",
                    "aptiva_acpafsiggfibrinogen_qual", "aptiva_acpafsigghistone1_qual", 
                    "aptivaccp3iga_qual", "aptiva_acpafsigavimentin2_qual",
                    "aptiva_acpafsigafibrinogen_qual", "aptiva_acpafsigahistone1_qual")

dat_2 <- dat_2_raw %>% 
  dplyr::rename(diagnosis=casecontrol, 
                study_id=masterstudyid,
                subj_id=masterstudyidnumeric,
                sampnum=sampordernum,
                age=ageserumsample,
                t_days=d_serum_ref) %>% 
  select(study_id, subj_id, sampnum,diagnosis, age, t_days, 
         year_dref, gender, eversmoke, familyhxra, race_ethnic, 
         all_of(biomarkers), all_of(bin_biomarkers))

dat_2$subj_id <- ifelse(dat_2$diagnosis=="Control", paste0(dat_2$subj_id, "_1"), paste(dat_2$subj_id)) # control and case had the same subj ID
dat_2$subj_id <- as.character(as.integer(factor(dat_2$subj_id)))
dat_2$study_id <- as.character(as.integer(factor(dat_2$study_id)))
dat_2$t_yrs <- dat_2$t_days/365 # changing from days to years

dat_2[,biomarkers] <- apply(dat_2[,biomarkers], 2, as.numeric)
dat_2 <- arrange(dat_2, as.numeric(subj_id), sampnum) %>% drop_na(all_of(biomarkers))

# Outcome matrix
Y <- as.matrix(subset(dat_2, select = biomarkers))
colnames(Y) <- c("aptivaccp3igg","aptiva_acpafsiggvimentin2",
                 "aptiva_acpafsiggfibrinogen","aptiva_acpafsigghistone1",
                 "aptivaccp3iga","aptiva_acpafsigavimentin2",
                 "aptiva_acpafsigafibrinogen","aptiva_acpafsigahistone1")
# Y <- apply(Y, 2, function(z) ifelse(z == 0, 1e-6, z))
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(dat_2$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
time <- dat_2$t_yrs # time before diagnosis
fem <- ifelse(dat_2$gender == "F", 1, 0) # indeicator for female
white <- ifelse(dat_2$race_ethnic == "W", 1, 0) # indicator for white
black <- ifelse(dat_2$race_ethnic == "B", 1, 0) # indicator for black
hispanic <- ifelse(dat_2$race_ethnic == "H", 1, 0) # indicator for hispanic
famhx <- ifelse(dat_2$familyhxra == "Yes", 1, 0)
bage.tmp <- dat_2 %>% 
  arrange(subj_id, t_yrs) %>%
  group_by(subj_id) %>% 
  filter(row_number() == 1) %>%
  mutate(bage = round(age - t_yrs)) %>%
  select(bage)
age_diag <- merge(dat_2, bage.tmp, by = "subj_id")$bage
diagnosis <- ifelse(dat_2$diagnosis == "Case", 1, 0)

rm(bage.tmp)

# ID vector
subj_id <- as.integer(factor(dat_2$subj_id, levels = unique(dat_2$subj_id)))
study_id <- as.integer(factor(dat_2$study_id, levels = unique(dat_2$study_id)))

maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)
D <- apply(Y, 2, function(z) as.numeric(z == max(z)))

clean <- data.frame(time, age_diag, fem, white, black, hispanic, famhx, subj_id, study_id, diagnosis, Y)

# Dichotomus biomarker values

Y_bin <- dat_2 %>%
  select(all_of(bin_biomarkers))

Y_bin <- as.matrix(Y_bin)
Y_bin <- apply(Y_bin, 2, as.integer)
colnames(Y_bin) <- colnames(Y)


