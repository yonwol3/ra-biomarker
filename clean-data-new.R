#-----------------------#
# dataset 2 cleaning
#-----------------------#
library(Microsoft365R)
library(readxl)
library(tidyverse)

# Kevin's loading code
#file_path <- "~/Dropbox/Projects/ra-biomarker/Data/KevinDat2.xlsx"
#dat_2raw <- read_xlsx(file_path)

# Yonatan's loading code
onedrive<- get_business_onedrive()
file_path <- "Attachments/KevinDat2.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
onedrive$download_file(
  src = file_path,
  dest = temp_file,
  overwrite = TRUE
)
dat_2raw<- read_xlsx(temp_file)
unlink(temp_file)
colnames(dat_2raw) <- tolower(colnames(dat_2raw))

biomarkers<-c("aptivaccp3iga_≥5#00flu","aptivaccp3igg_≥5#00flu",
              "aptivapad1igg_≥5#00au","aptivapad4igg_≥5#00au",
              "aptiva_acpafsiggvimentin2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
              "aptivapad1iga_≥5#00au","aptivapad4iga_≥5#00au",
              "aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")

dat_2 <- dat_2raw %>% 
  dplyr::rename(diagnosis=casecontrol, 
                study_id=masterstudyid,
                subj_id=masterstudyidnumeric,
                sampnum=sampordernum,
                age=ageserumsample,
                t_days=d_serum_ref) %>% 
  select(study_id, subj_id, sampnum,diagnosis, age,t_days,year_dref,
         gender,eversmoke,familyhxra,race_ethnic, all_of(biomarkers))

dat_2$subj_id <- ifelse(dat_2$diagnosis=="Control", paste0(dat_2$subj_id, "_1"), paste(dat_2$subj_id)) # control and case had the same subj ID
dat_2$subj_id <- as.character(as.integer(factor(dat_2$subj_id)))
dat_2$t_yrs <- dat_2$t_days/365 # changing from days to years

dat_2[,biomarkers] <- apply(dat_2[,biomarkers], 2, as.numeric)
dat_2 <- arrange(dat_2, as.numeric(subj_id), sampnum) %>% drop_na(all_of(biomarkers))

# Outcome matrix
Y <- as.matrix(subset(dat_2, select = biomarkers))
Y <- apply(Y, 2, function(z) ifelse(z == 0, 1e-6, z))
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(dat_2$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

# Covariates
time <- dat_2$t_yrs # time before diagnosis
fem <- ifelse(dat_2$gender == "F", 1, 0) # indeicator for female
nw <- ifelse(dat_2$race_ethnic == "W", 0, 1) # indicator for non-white
famhx <- ifelse(dat_2$familyhxra == "No", 0, 1)
bage.tmp <- dat_2$age[unique(dat_2$subj_id)]
bage <- rep(bage.tmp, times = table(dat_2$subj_id))
diagnosis <- ifelse(dat_2$diagnosis == "RA", 1, 0)

rm(bage.tmp)

# ID vector
subj_id <- as.integer(factor(dat_2$subj_id, levels = unique(dat_2$subj_id)))
study_id <- as.integer(factor(dat_2$study_id, levels = unique(dat_2$study_id)))

cens_max <- apply(logY, 2, function(z) as.numeric(z == max(z, na.rm = T))) # binary yes/no
cens_min <- apply(logY, 2, function(z) as.numeric(z == min(z, na.rm = T))) # binary yes/no
maxY <- apply(logY, 2, max, na.rm = T) # max value for each biomarker vector
minY <- apply(logY, 2, min, na.rm = T) # min value for each biomarker vector

#----------------------------------------------------------------#
#Fun that outputs stan objects from input of a pair of biomarker
#----------------------------------------------------------------#

create_biomarker_objects <- function(iga_biomarker, igg_biomarker, data = dat_2) {
  
  strip_suffix <- function(x) {
    x_no_suffix <- sub("_≥5#00.*", "", x)  # Remove everything after '_≥5#00'
    x_no_suffix <- sub("iga", "", x_no_suffix) # Remove 'iga'
    x_no_suffix <- sub("igg", "", x_no_suffix) # Remove 'igg'
    x_no_suffix
  }
  
  XX_iga <- strip_suffix(iga_biomarker)
  XX_igg <- strip_suffix(igg_biomarker)
  
  if (XX_iga != XX_igg) {
    stop("The IgA and IgG biomarker names do not share the same root identifier.")
  }
  XX <- XX_iga
  
  if (!all(c(iga_biomarker, igg_biomarker) %in% names(data))) {
    stop("The provided biomarker columns do not exist in the dataset.")
  }
  
  Y_XX <- as.matrix(data[, c(iga_biomarker, igg_biomarker)])
  Y_XX <- apply(Y_XX, 2, function(z) ifelse(z == 0, 1e-6, z))
  
  logY_XX <- log(Y_XX)
  N_XX <- nrow(Y_XX)                       
  M_XX <- nlevels(factor(data$subj_id))    
  K_XX <- ncol(Y_XX)                       
  
  cens_max_XX <- apply(logY_XX, 2, function(z) as.numeric(z == max(z, na.rm = TRUE)))
  cens_min_XX <- apply(logY_XX, 2, function(z) as.numeric(z == min(z, na.rm = TRUE)))
  
  maxY_XX <- apply(logY_XX, 2, max, na.rm = TRUE)
  minY_XX <- apply(logY_XX, 2, min, na.rm = TRUE)
  
  
  result_list <- list(
    Y = Y_XX,
    logY = logY_XX,
    N = N_XX,
    M = M_XX,
    K = K_XX,
    cens_max = cens_max_XX,
    cens_min = cens_min_XX,
    maxY = maxY_XX,
    minY = minY_XX
  )
  
  # Name the list as XX
  names(result_list) <- paste0(names(result_list), "_", XX)
  
  return(result_list)
}

