#-----------------------#
# dataset 2 cleaning
#-----------------------#

library(Microsoft365R)
library(readxl)
library(tidyverse)
onedrive<- get_business_onedrive()
file_path <- "Attachments/KevinDat2.xlsx"

onedrive$download_file(
  src = file_path,
  dest = temp_file,
  overwrite = TRUE
)

dat_2raw<- read_xlsx(temp_file)
unlink(temp_file)

colnames(dat_2raw)<-tolower(colnames(dat_2raw))


biomarkers<-c("aptivaccp3iga_≥5#00flu","aptivaccp3igg_≥5#00flu","aptivapad1igg_≥5#00au","aptivapad2igg_≥5#00au",
              "aptivapad3igg_≥5#00au","aptivapad4igg_≥5#00au","aptivapad6igg_≥5#00au","aptiva_acpafsiggvimentin2_≥5#00au",
              "aptiva_acpafsigghistone2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au","aptiva_acpafsiggvimentin1_≥5#00au",
              "aptivapad1iga_≥5#00au","aptivapad2iga_≥5#00au","aptivapad3iga_≥5#00au","aptivapad4iga_≥5#00au",
              "aptivapad6iga_≥5#00au","aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigahistone2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au",
              "aptiva_acpafsigahistone1_≥5#00au","aptiva_acpafsigavimentin1_≥5#00au","quantalitecarpigg_≥20au",
              "quantalitecarpiga_≥20au","quantaflashrfiga_≥20#0cu_a","quantaflashrfigm_≥5#0iuml_a")
dat_2<-dat_2raw %>% 
        dplyr::rename(diagnosis=casecontrol, 
               study_id=masterstudyid,
               subj_id=masterstudyidnumeric,
               sampnum=sampordernum,
               age=ageserumsample,
               t_days=d_serum_ref) %>% 
        select(study_id, subj_id, sampnum,diagnosis, age,t_days,year_dref,
               gender,eversmoke,familyhxra,race_ethnic,all_of(biomarkers)
               ) 
dat_2$subj_id<-ifelse(dat_2$diagnosis=="Control", paste0(dat_2$subj_id,"_1"),paste(dat_2$subj_id))# control and case had the same subj ID
dat_2$subj_id <- as.character(as.integer(factor(dat_2$subj_id)))
dat_2$t_yrs<-dat_2$t_days/365 # changing from days to years

mean_age<-dat_2 %>% arrange(subj_id,sampnum) %>% select(subj_id,age) 
mean_age$subj_id<-as.character(mean_age$subj_id)
mean_age<- mean_age %>% 
             dplyr::group_by(subj_id) %>% 
             dplyr::summarise(mean_age=mean(age))
dat_2<-dplyr::left_join(dat_2,mean_age, by="subj_id") #mean age created
dat_2<-arrange(dat_2,as.numeric(subj_id),sampnum)

# Outcome matrix
Y <- as.matrix( subset(dat_2, select = biomarkers) )
Y[] <- lapply(Y, as.numeric)
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
subj_id <- dat_2$subj_id
study_id <- dat_2$study_id

cens_max <- apply(logY, 2, function(z) as.numeric(z == max(z)))
cens_min <- apply(logY, 2, function(z) as.numeric(z == min(z)))
maxY <- apply(logY, 2, max)
minY <- apply(logY, 2, min)




  