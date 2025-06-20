#############################################
## PURPOSE: Script for creation of tables  ##
##          in RA data analysis            ##
## BY:      Yonatan Wolde                 ##
#############################################

# Libraries


library(tidyverse)
library(labelled)
library(table1)
library(naniar)
library(mcmcse)
library(gt)
library(kableExtra)

## Table 1

# Load the dataframes
source("clean-data.R")
source("clean-data-new.R")

preserve_factor <- function(x) {
                      lab <- attr(x, "label")
                      fx <- factor(x)
                      attr(fx, "label") <- lab
                      fx
                    }


# Add age at diagnosis values to each dataset
clean$agediag <- raDat$age - raDat$t_yrs # adding age at diagnosis from raDat
dat_2$agediag <- dat_2$age - dat_2$t_yrs  # age at diagnosis in new dataset

# Rename covariate elements for the 'clean' dataset
clean$fem <- factor(clean$fem, labels = c("Male", "Female"))
clean$nw <- factor(clean$white, labels = c("Non-white","White"))
clean$famhx<-raDat$familyhxra 
clean$eversmoke <- raDat$eversmoke  # smoking status from raDat
clean$diagnosis <- factor(clean$diagnosis, labels = c("No RA", "RA"))
clean <- clean %>%
  dplyr::mutate(
    famhx_c = dplyr::case_when(
      famhx == "Yes" ~ "Yes",
      famhx == "No" ~ "No",
      T ~ "Not available"
    ),
    smoking = dplyr::case_when(
      eversmoke == "No" ~ "No",
      eversmoke == "Yes" ~ "Yes",
      TRUE ~ "Not available"
    )
  )

# Rename covariate elements for the 'dat_2' dataset
dat_2$gender <- factor(dat_2$gender, levels = c("M", "F"), labels = c("Male", "Female"))
dat_2$nw <- ifelse(dat_2$race_ethnic == "W", 0, 1)  # non-white indicator
dat_2$nw <- factor(dat_2$nw, labels = c("White", "Non-white"))
dat_2$diagnosis <- factor(dat_2$diagnosis, 
                          levels = c("Control", "Case"), 
                          labels = c("No RA", "RA"))
dat_2$subj_id<-subj_id
dat_2$study_id<-study_id
dat_2 <- dat_2 %>%
  dplyr::mutate(
    famhx_c = dplyr::case_when(
      familyhxra == "No" ~ "No",
      familyhxra == "Yes" ~ "Yes",
      TRUE ~ "Not available"
    ),
    smoking = dplyr::case_when(
      eversmoke == "No" ~ "No",
      eversmoke == "Yes" ~ "Yes",
      TRUE ~ "Not available"
    )
  )

# Compute mean biomarker values per subject to create Table 1
clean_tb1 <- clean %>% 
  dplyr::group_by(subj_id) %>%
  dplyr::summarise(
    gender = dplyr::first(fem),
    race = dplyr::first(nw),
    famhx_c = dplyr::first(famhx_c),
    diagnosis = dplyr::first(diagnosis),
    smoking = dplyr::first(smoking),
    mean_agediag = mean(agediag),
    dplyr::across(starts_with("ig"), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}")
  )

dat_2_tb1 <- dat_2 %>% 
  dplyr::group_by(subj_id) %>%
  dplyr::summarise(
    gender = dplyr::first(gender),
    race = dplyr::first(nw),
    famhx_c = dplyr::first(famhx_c),
    smoking = dplyr::first(smoking),
    diagnosis = dplyr::first(diagnosis),
    mean_agediag = mean(agediag),
    dplyr::across(starts_with("apti"), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}")
  )

# Specify variable labels for the final Table 1 output
clean_tb1 <- set_variable_labels(clean_tb1,
                                 gender = "Sex",
                                 race = "Race",
                                 famhx_c = "Family RA History",
                                 smoking = "Eversmoke", 
                                 mean_agediag = "Age at Diagnosis",
                                 mean_igarfconc_ = "RF IgA",
                                 mean_igmrfconc_ = "RF IgM",
                                 mean_iggrfconc_ = "RF IgG",
                                 mean_igaccpavgconc = "ACPA IgA",
                                 mean_igmccpavgconc = "ACPA IgM",
                                 mean_iggccpavgconc = "ACPA IgG")





dat_2_tb1 <- set_variable_labels(dat_2_tb1,
                                 gender = "Sex",
                                 race = "Race",
                                 famhx_c = "Family RA History",
                                 smoking = "Eversmoke", 
                                 mean_agediag = "Age at Diagnosis",
                                 `mean_aptivaccp3iga_≥5#00flu`="ccp3iga",
                                 `mean_aptivaccp3igg_≥5#00flu`="ccp3igg",
                                 `mean_aptivapad1igg_≥5#00au`="pad1igg",             
                                 `mean_aptivapad4igg_≥5#00au`="pad4igg",
                                 `mean_aptiva_acpafsiggvimentin2_≥5#00au`= "acpafsiggvimentin2",
                                 `mean_aptiva_acpafsiggfibrinogen_≥5#00au`="acpafsiggfibrinogen",
                                 `mean_aptiva_acpafsigghistone1_≥5#00au`= "acpafsigghistone1" ,
                                 `mean_aptivapad1iga_≥5#00au`="pad1iga", 
                                 `mean_aptivapad4iga_≥5#00au`="pad4iga",           
                                 `mean_aptiva_acpafsigavimentin2_≥5#00au`="acpafsigavimentin2" , 
                                 `mean_aptiva_acpafsigafibrinogen_≥5#00au`="acpafsigafibrinogen",
                                 `mean_aptiva_acpafsigahistone1_≥5#00au`="acpafsigahistone1" )






#--------------------#
#
# table-1 BOTH cohorts
#
#---------------------#
# Set variable labels for dat_2_tb1

# Create Table 1 datasets by removing the subject ID column and renaming 'diagnosis' as 'group'
df1 <- clean_tb1[, -1] %>% 
  dplyr::rename(group = diagnosis)
df2 <- dat_2_tb1[, -1] %>% 
  dplyr::rename(group = diagnosis)


# Make sure group is a factor with proper labels.
df1$group <- factor(df1$group, labels = c("Non-RA", "RA"))
df2$group <- factor(df2$group, labels = c("Non-RA", "RA"))

# Helper: get the label of a variable (if it exists) or return the variable name.
get_label <- function(data, var) {
  lab <- attr(data[[var]], "label")
  if (is.null(lab)) var else lab
}

# Summaries for numeric variables: Continuous variables summarized as Median (IQR)
summarize_numeric_vars <- function(data, numeric_vars, group_var) {
  out_list <- lapply(numeric_vars, function(v) {
    # Calculate summary stats by group.
    stat_df <- data %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::summarize(
        med_val = median(.data[[v]], na.rm = TRUE),
        q1 = quantile(.data[[v]], probs = 0.25, na.rm = TRUE),
        q3 = quantile(.data[[v]], probs = 0.75, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        MedIQR = sprintf("%.2f (%.2f-%.2f)", med_val, q1, q3)
      ) %>%
      dplyr::select(dplyr::all_of(group_var), MedIQR)
    
    # Look up the label for the variable.
    label_name <- get_label(data, v)
    
    stat_long <- stat_df %>%
      dplyr::mutate(stat = "Median (IQR)", var = label_name)
    
    # Rename the computed value column to "val"
    names(stat_long)[names(stat_long) == "MedIQR"] <- "val"
    stat_long
  })
  
  all_numeric <- dplyr::bind_rows(out_list) %>%
    tidyr::pivot_wider(
      id_cols = c("var", "stat"),
      names_from = group_var,
      values_from = val
    )
  
  return(all_numeric)
}

# Summaries for categorical variables
summarize_categorical_vars <- function(data, factor_vars, group_var) {
  out_list <- lapply(factor_vars, function(v) {
    dcat <- data %>%
      dplyr::group_by(.data[[group_var]], .data[[v]]) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::mutate(percent = round(100 * n / sum(n), 1)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(val = paste0(n, " (", percent, "%)")) %>%
      dplyr::select(dplyr::all_of(group_var), level = dplyr::all_of(v), val)
    
    # Look up the label for the variable.
    label_name <- get_label(data, v)
    
    dcat_long <- dcat %>%
      dplyr::rename(stat = level) %>%
      dplyr::mutate(var = label_name)
    
    dcat_wide <- dcat_long %>%
      tidyr::pivot_wider(
        id_cols = c("var", "stat"),
        names_from = group_var,
        values_from = "val"
      )
    
    dcat_wide
  })
  
  all_cats <- dplyr::bind_rows(out_list)
  return(all_cats)
}

# One-sample summarizer: generates the summary table for one dataset.
make_summary_one_sample <- function(data, group_var = "group") {
  numeric_vars <- names(data)[sapply(data, is.numeric)]
  factor_vars <- names(data)[sapply(data, function(x) is.character(x) | is.factor(x))]
  
  numeric_vars <- setdiff(numeric_vars, group_var)
  factor_vars <- setdiff(factor_vars, group_var)
  
  # Ensure factor variables are factors.
  data[factor_vars] <- lapply(data[factor_vars], preserve_factor)
  
  df_numeric <- NULL
  if (length(numeric_vars) > 0) {
    df_numeric <- summarize_numeric_vars(data, numeric_vars, group_var)
  }
  
  df_cat <- NULL
  if (length(factor_vars) > 0) {
    df_cat <- summarize_categorical_vars(data, factor_vars, group_var)
  }
  
  combined <- dplyr::bind_rows(df_cat, df_numeric)
  return(combined)
}

# Create summaries for your two samples.
summary_df1 <- make_summary_one_sample(df1, "group") %>%
  dplyr::rename(Sample1_RA = RA, Sample1_NonRA = `Non-RA`)

summary_df2 <- make_summary_one_sample(df2, "group") %>%
  dplyr::rename(Sample2_RA = RA, Sample2_NonRA = `Non-RA`)

# Combine the two summaries; replace NA with "-"
combined_table <- dplyr::full_join(summary_df1, summary_df2, by = c("var", "stat")) %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.), "-", .)))

# Compute RA / Non-RA sample sizes for each dataset.
n_sample1_RA <- sum(df1$group == "RA")
n_sample1_NonRA <- sum(df1$group == "Non-RA")
n_sample2_RA <- sum(df2$group == "RA")
n_sample2_NonRA <- sum(df2$group == "Non-RA")

# Build the final GT Table with bold formatting.
final_gt_table <- combined_table %>%
  gt::gt(
    rowname_col = "stat",    # row labels (e.g., "Median (IQR)", "Female", etc.)
    groupname_col = "var"      # groups each variable together (displayed in bold)
  ) %>%
  gt::tab_spanner(
    label = "Sample 1",
    columns = c("Sample1_RA", "Sample1_NonRA")
  ) %>%
  gt::tab_spanner(
    label = "Sample 2",
    columns = c("Sample2_RA", "Sample2_NonRA")
  ) %>%
  # Bold the RA/Non-RA column labels with sample sizes.
  gt::cols_label(
    Sample1_RA = gt::md(sprintf("**RA (N=%d)**", n_sample1_RA)),
    Sample1_NonRA = gt::md(sprintf("**Non-RA (N=%d)**", n_sample1_NonRA)),
    Sample2_RA = gt::md(sprintf("**RA (N=%d)**", n_sample2_RA)),
    Sample2_NonRA = gt::md(sprintf("**Non-RA (N=%d)**", n_sample2_NonRA))
  ) %>%
  gt::tab_header(
    title = gt::md("**Table 1: Descriptive Summary of Sample 1 and Sample 2**")
  ) %>%
  # Bold the group labels (variable names).
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups(groups = gt::everything())
  ) %>%
  # Add lines and center columns.
  gt::opt_table_lines(extent = "all") %>%
  gt::cols_align("center")






# --- Generate the final table using kable and kableExtra ---
var_labels<-unique(combined_table$var)
table1_obj <- combined_table %>%
  select(-var) %>%
  kable(
    caption = "Table 1: Descriptive Summary of Sample A and Sample B",
    col.names = rep("", ncol(.)),
    format = "pipe"
  ) %>%
  add_header_above(c(" " = 1, "Non-RA (N=210)" = 1, "RA (N=214)" = 1, "Non-RA (N=309)" = 1, "RA (N=309)" = 1)) %>%
  add_header_above(c(" " = 1, "Sample A" = 2, "Sample B" = 2)) %>%
  kable_styling(full_width = FALSE, position = "center")

for (i in 1:length(var_labels)) {
  if (i == 1) {
    table1_obj <- pack_rows(table1_obj, group_label = var_labels[i], start_row = 1, end_row = 2)
  } else if (i == 2) {
    table1_obj <- pack_rows(table1_obj, group_label = var_labels[i], start_row = 3, end_row = 4)
  } else if (i == 3) {
    table1_obj <- pack_rows(table1_obj, group_label = var_labels[i], start_row = 5, end_row = 7)
  } else if (i == 4) {
    table1_obj <- pack_rows(table1_obj, group_label = var_labels[i], start_row = 8, end_row = 10)
  } else {
    table1_obj <- pack_rows(table1_obj, group_label = var_labels[i], start_row = i + 6, end_row = i + 6)
  }
}
saveRDS(combined_table, file = "../../table1_obj_df.rds")
saveRDS(table1_obj, file = "../../table1_obj.rds")






#------------------------------------------------#
#
# change-point+ magnitude(90% prob being >0)
# Original Biomarkers
#
#------------------------------------------------#

# function to find highest posterior density interval

# hpd <- function(x, alpha = 0.05){
#   
#   # x is the vector of bootstrap estimates
#   n = length(x)
#   m = round(n * alpha)
#   x = sort(x)
#   y = x[(n - m + 1):n] - x[1:m]
#   z = min(y)
#   k = which(y == z)[1]
#   c(x[k], x[n - m + k])
#   
# }
# 
# 
# onedrive<- get_business_onedrive()
# file_path <- "Attachments/mcmc_trunc.RData"
# temp_file <- tempfile(fileext = ".RData")
# onedrive$download_file(
#   src = file_path,
#   dest = temp_file,
#   overwrite = TRUE
# )
# mcmc_trunc<- get(load(temp_file)[1])
# mcmc<-mcmc_trunc[[1]]
# 
# # Number of biomarkers being analyzed
# k <- 6
# 
# # Define the time grid
# time_grid <- seq(-20, 10, by = 0.01)
# 
# # Extract kappa and gamma from the mcmc data
# kappa <- mcmc[, 7:12]
# gamma <- mcmc[, 1:6]
# 
# biomarker_labels <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")
# time_labels <- as.character(time_grid)
# iteration_labels <- paste0("iter", seq_len(nrow(kappa)))
# 
# # Create a 3-dimensional array with named dimensions:
# # Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
# res <- array(NA, dim = c(k, length(time_grid), nrow(kappa)),
#              dimnames = list(
#                biomarker = biomarker_labels,
#                time = time_labels,
#                iteration = iteration_labels
#              ))
# 
# # Loop over biomarkers, time grid, and iterations to fill the array.
# for (b in 1:k) {
#   gamma_b <- gamma[, b]     # gamma values for biomarker b
#   kappa_b <- kappa[, b]     # kappa values for biomarker b
#   
#   for (i in seq_along(time_grid)) {
#     t <- time_grid[i]
#     for (j in 1:nrow(kappa)) {
#       res[b, i, j] <- (t - kappa_b[j]) * gamma_b[j]
#     }
#   }
# }
# 
# # Convert the array to a data frame.
# result_df <- as.data.frame.table(res, responseName = "value")
# colnames(result_df) <- c("biomarker", "time", "iteration", "value")
# 
# result_df$time <-as.numeric(as.character(result_df$time))
# 
# # Calculate the proportion of iterations for which 'value' > 0, for each biomarker and time point.
# result_df <- result_df %>% 
#   dplyr::group_by(biomarker, time) %>% 
#   dplyr::summarise(prop_positive = mean(value > 0),.groups = "drop") %>% 
#   arrange(desc(prop_positive))
# 
# # select time where absolute difference between prop_positive and 0.9 is minimized
# 
# closest_threshold <- result_df %>% 
#   dplyr::group_by(biomarker) %>% 
#   dplyr::slice(which.min(abs(prop_positive - 0.9))) %>% 
#   dplyr::ungroup()%>% 
#   arrange(time)
# 
# 
# # Gamma only (mean (95% HPDI))
# 
# gamma_summ<-matrix(NA,ncol = 2, nrow=length(biomarker_labels))
# for ( i in 1:length(biomarker_labels)) {
#   gamma_b<-gamma[, i]
#   mean<-round(mean(gamma_b),2)
#   q_2<- round(hpd(gamma_b)[1], 2) # 2.5th quantile
#   q_97<-round(hpd(gamma_b)[2], 2) # 97.5th quantile
#   gamma_summ[i, 2]<-paste(mean, "[",q_2,"-",q_97,"]")
#   gamma_summ[i, 1]<-biomarker_labels[i]
# }
# gamma_summ<- as.data.frame(gamma_summ)
# colnames(gamma_summ)<- c("biomarker","gamma mean[95% CI]")
# 
# # kappa only (mean (95% HPDI))
# kappa_summ <- data.frame(biomarker=character(0),kappa_summ=numeric(0))
# for (i in 1:ncol(kappa)) {
#   kappa_b<-kappa[, i]
#   mean<-round(mean(kappa_b),2)
#   q_2<- round(hpd(kappa_b)[1], 2) # 2.5th quantile
#   q_97<-round(hpd(kappa_b)[2], 2) 
#   kappa_summ[i, 2]<-paste(mean, "[",q_2,"-",q_97,"]")
#   kappa_summ[i, 1]<-biomarker_labels[i]
# }
# 
# colnames(kappa_summ)<- c("biomarker","kappa mean[95% CI]")
# 
# closest_threshold<-left_join(closest_threshold,kappa_summ, by="biomarker") 
# closest_threshold<-left_join(closest_threshold, gamma_summ, by="biomarker")
# write.csv(closest_threshold,"../../original_prob_ci.csv")
# 
# #------------------------------------------------#
# #
# # change-point+ magnitude(90% prob being >0)
# # New Biomarkers
# #
# #------------------------------------------------#
# 
# onedrive<- get_business_onedrive()
# file_path <- "Attachments/mcmc_new_trunc.RData"
# temp_file <- tempfile(fileext = ".RData")
# onedrive$download_file(
#   src = file_path,
#   dest = temp_file,
#   overwrite = TRUE
# )
# mcmc_trunc_new<- get(load(temp_file)[1])
# mcmc_new<-mcmc_trunc_new[[1]]
# 
# # Number of biomarkers being analyzed
# k <- 12
# 
# # Define the time grid
# time_grid <- seq(-20, 10, by = 0.01)
# 
# # Extract kappa and gamma from the mcmc data
# kappa <- mcmc_new[, 13:24]
# gamma <- mcmc_new[, 1:12]
# 
# # Define labels for each dimension
# biomarkers<-c("aptivaccp3iga","aptivaccp3igg",
#               "aptivapad1igg_≥5#00au","aptivapad4igg_≥5#00au",
#               "aptiva_acpafsiggvimentin2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
#               "aptivapad1iga_≥5#00au","aptivapad4iga_≥5#00au",
#               "aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
# biomarkers <- sub("^aptiva", "", biomarkers)
# biomarkers<- sub("_≥5#00au", "", biomarkers)
# biomarkers<-sub("_","",biomarkers)
# biomarker_labels <- biomarkers
# 
# time_labels <- as.character(time_grid)
# iteration_labels <- paste0("iter", seq_len(nrow(kappa)))
# 
# # Create a 3-dimensional array with named dimensions:
# # Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
# res <- array(NA, dim = c(k, length(time_grid), nrow(kappa)),
#              dimnames = list(
#                biomarker = biomarker_labels,
#                time = time_labels,
#                iteration = iteration_labels
#              ))
# 
# # Loop over biomarkers, time grid, and iterations to fill the array.
# for (b in 1:k) {
#   gamma_b <- gamma[, b]     # gamma values for biomarker b
#   kappa_b <- kappa[, b]     # kappa values for biomarker b
#   
#   for (i in seq_along(time_grid)) {
#     t <- time_grid[i]
#     for (j in 1:nrow(kappa)) {
#       res[b, i, j] <- (t - kappa_b[j]) * gamma_b[j]
#     }
#   }
# }
# 
# # Convert the array to a data frame.
# result_df_new<- as.data.frame.table(res, responseName = "value")
# colnames(result_df_new) <- c("biomarker", "time", "iteration", "value")
# 
# result_df_new$time <-as.numeric(as.character(result_df_new$time))
# 
# # Calculate the proportion of iterations for which 'value' > 0, for each biomarker and time point.
# result_df_new <- result_df_new %>% 
#   dplyr::group_by(biomarker, time) %>% 
#   dplyr::summarise(prop_positive = mean(value > 0),.groups = "drop") %>% 
#   arrange(desc(prop_positive))
# 
# # select time where absolute difference between prop_positive and 0.9 is minimized
# 
# closest_threshold_new<- result_df_new %>% 
#   dplyr::group_by(biomarker) %>% 
#   dplyr::slice(which.min(abs(prop_positive - 0.9))) %>% 
#   dplyr::ungroup()%>% 
#   arrange(time)
# 
# # Gamma only (mean (95% HPDI))
# 
# gamma_summ<-matrix(NA,ncol = 2, nrow=length(biomarker_labels))
# for ( i in 1:length(biomarker_labels)) {
#   gamma_b<-gamma[, i]
#   mean<-round(mean(gamma_b),2)
#   q_2<- round(hpd(gamma_b)[1], 2) # 2.5th quantile
#   q_97<-round(hpd(gamma_b)[2], 2) # 97.5th quantile
#   gamma_summ[i, 2]<-paste(mean, "[",q_2,"-",q_97,"]")
#   gamma_summ[i, 1]<-biomarker_labels[i]
# }
# gamma_summ<- as.data.frame(gamma_summ)
# colnames(gamma_summ)<- c("biomarker","gamma mean[95% CI]")
# 
# # kappa only (mean (95% HPDI))
# kappa_summ <- data.frame(biomarker=character(0),kappa_summ=numeric(0))
# for (i in 1:ncol(kappa)) {
#   kappa_b<-kappa[, i]
#   mean<-round(mean(kappa_b),2)
#   q_2<- round(hpd(kappa_b)[1], 2) # 2.5th quantile
#   q_97<-round(hpd(kappa_b)[2], 2) 
#   kappa_summ[i, 2]<-paste(mean, "[",q_2,"-",q_97,"]")
#   kappa_summ[i, 1]<-biomarker_labels[i]
# }
# 
# colnames(kappa_summ)<- c("biomarker","kappa mean[95% CI]")
# 
# closest_threshold_new<-left_join(closest_threshold_new,kappa_summ, by="biomarker") 
# closest_threshold_new<-left_join(closest_threshold_new, gamma_summ, by="biomarker")
# write.csv(closest_threshold_new,"../../new_prob_ci.csv")


