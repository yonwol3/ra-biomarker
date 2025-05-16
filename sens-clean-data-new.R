source ("clean-data-new.R")

# remove observations that are t>0 (after diagnosis)
# remove from both cases and controls

dat2_sens<- dat_2 %>% 
  filter(t_yrs<=0)

Y <- as.matrix(dat2_sens[, 12:23])
logY <- log(Y) # log transform responses

# Sample numbers
N <- nrow(Y) # number of samples
M <- nlevels(factor(dat2_sens$subj_id)) # number of participants
K <- ncol(Y) # number of measurements per sample

cens_max <- apply(Y, 2, function(z) as.numeric(z == max(z)))
cens_min <- apply(Y, 2, function(z) as.numeric(z == min(z)))
maxY <- apply(Y, 2, max)
minY <- apply(Y, 2, min)
L <- matrix(rep(minY, N), ncol = K, byrow = TRUE)
U <- matrix(rep(maxY, N), ncol = K, byrow = TRUE)

# Dichotomous outcome dataset 

onedrive<- get_business_onedrive()
file_path <- "Attachments/KevinDat2.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
onedrive$download_file(
  src = file_path,
  dest = temp_file,
  overwrite = TRUE
)

dat_2_diraw <- read_xlsx(temp_file)
unlink(temp_file)
colnames(dat_2_diraw) <- tolower(colnames(dat_2_diraw))

biomarkers<-c("aptivaccp3iga_≥5#00flu","aptivaccp3igg_≥5#00flu",
              "aptivapad1igg_≥5#00au","aptivapad4igg_≥5#00au",
              "aptiva_acpafsiggvimentin2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
              "aptivapad1iga_≥5#00au","aptivapad4iga_≥5#00au",
              "aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
idx<-which(colnames(dat_2_diraw) %in% biomarkers) 
idx<- idx +1

# Dichotomus biomarker values

di_biomarkers<-colnames(dat_2_diraw)[idx]
a<- dat_2_diraw %>%
  select(biomarkers) %>% 
  drop_na()
dat_2_di <- dat_2_diraw %>% 
  dplyr::rename(diagnosis=casecontrol, 
                study_id=masterstudyid,
                subj_id=masterstudyidnumeric,
                sampnum=sampordernum,
                age=ageserumsample,
                t_days=d_serum_ref) %>% 
  select(study_id, subj_id, sampnum,diagnosis, age,t_days,year_dref,
         gender,eversmoke,familyhxra,race_ethnic, all_of(di_biomarkers))

dat_2_di$subj_id <- ifelse(dat_2_di$diagnosis=="Control", paste0(dat_2_di$subj_id, "_1"), paste(dat_2_di$subj_id)) # control and case had the same subj ID
dat_2_di$subj_id <- as.character(as.integer(factor(dat_2_di$subj_id)))
dat_2_di$t_yrs <- dat_2_di$t_days/365 # changing from days to years
dat_2_di[ ,di_biomarkers] <- apply(dat_2_di[,di_biomarkers], 2, as.numeric)
dat_2_di <- arrange(dat_2_di, as.numeric(subj_id), sampnum) %>% drop_na(all_of(di_biomarkers))



