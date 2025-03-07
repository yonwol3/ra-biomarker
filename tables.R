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

## Table 1


# Load the dataframes
source("clean-data.R")
source("clean-data-new.R")

# Add age at diagnosis values to each dataset
clean$agediag <- raDat$agediag  # adding mean age at diagnosis from raDat
dat_2$agediag <- dat_2$age - dat_2$t_yrs  # age at diagnosis in new dataset

# Rename covariate elements for the 'clean' dataset
clean$fem <- factor(clean$fem, labels = c("Male", "Female"))
clean$nw <- factor(clean$nw, labels = c("White", "Non-white"))
clean$eversmoke <- raDat$eversmoke  # smoking status from raDat
clean$diagnosis <- factor(clean$diagnosis, labels = c("No RA", "RA"))
clean <- clean %>%
  dplyr::mutate(
    famhx_c = dplyr::case_when(
      famhx == 0 ~ "No",
      famhx == 1 ~ "Yes",
      is.na(famhx) ~ "Not available"
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

##### Table-1 just for the original biomarkers for now

render.median.IQR <- function(x, ...) {
  c('', 
    `Median (IQR)` = sprintf("%s (%s, %s)", round(median(x),1), 
                             round(quantile(x, 0.25),1), round(quantile(x, 0.75),1)))
}
# Define the variables:
covariates <- c("gender", "race", "famhx_c", "smoking")
biomarkers <- c("mean_agediag", "mean_igarfconc_", "mean_igmrfconc_", 
                "mean_iggrfconc_", "mean_igaccpavgconc", "mean_igmccpavgconc", 
                "mean_iggccpavgconc")
t1_variables <- c(covariates, biomarkers)

# If needed, enclose variable names with backticks:
t1_variables_backticked <- paste0("`", t1_variables, "`")

# Create the table with table1 using our custom continuous renderer:
t1 <- table1(
  as.formula(paste0("~", paste(t1_variables_backticked, collapse = "+"), "| diagnosis")),
  data = clean_tb1,
  caption = "Table1: Descriptive summary of original sample by RA status",
  render.cont = render.median.IQR
)

t1



#--------------------#
#
# table-1 BOTH cohorts
#
#---------------------#

# Create Table 1 datasets by removing the subject ID column and renaming 'diagnosis' as 'group'
df1 <- clean_tb1[, -1] %>% 
  dplyr::rename(group = diagnosis)
df2 <- dat_2_tb1[, -1] %>% 
  dplyr::rename(group = diagnosis)
df1$group <- factor(df1$group, labels = c("Non-RA", "RA"))
df2$group <- factor(df2$group, labels = c("Non-RA", "RA"))

# Summaries for numeric variables: Continuous variables summarized as Median (IQR)
summarize_numeric_vars <- function(data, numeric_vars, group_var) {
  out_list <- lapply(numeric_vars, function(v) {
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
    
    stat_long <- stat_df %>%
      dplyr::mutate(stat = "Median (IQR)", var = v)
    # Use base R to rename 'MedIQR' to 'val'
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
    
    dcat_long <- dcat %>%
      dplyr::rename(stat = level) %>%
      dplyr::mutate(var = v)
    
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

# One-sample summarizer
make_summary_one_sample <- function(data, group_var = "group") {
  numeric_vars <- names(data)[sapply(data, is.numeric)]
  factor_vars <- names(data)[sapply(data, function(x) is.character(x) | is.factor(x))]
  numeric_vars <- setdiff(numeric_vars, group_var)
  factor_vars <- setdiff(factor_vars, group_var)
  
  data[factor_vars] <- lapply(data[factor_vars], factor)
  
  df_numeric <- NULL
  if (length(numeric_vars) > 0) {
    df_numeric <- summarize_numeric_vars(data, numeric_vars, group_var)
  }
  
  df_cat <- NULL
  if (length(factor_vars) > 0) {
    df_cat <- summarize_categorical_vars(data, factor_vars, group_var)
  }
  
  combined <- dplyr::bind_rows(df_numeric, df_cat)
  return(combined)
}

# Create summaries for your two samples
summary_df1 <- make_summary_one_sample(df1, "group") %>%
  dplyr::rename(Sample1_RA = RA, Sample1_NonRA = `Non-RA`)

summary_df2 <- make_summary_one_sample(df2, "group") %>%
  dplyr::rename(Sample2_RA = RA, Sample2_NonRA = `Non-RA`)

combined_table <- dplyr::full_join(summary_df1, summary_df2, by = c("var", "stat")) %>%
  dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.), "-", .)))

# Compute RA / Non-RA sample sizes for each dataset
n_sample1_RA <- sum(df1$group == "RA")
n_sample1_NonRA <- sum(df1$group == "Non-RA")
n_sample2_RA <- sum(df2$group == "RA")
n_sample2_NonRA <- sum(df2$group == "Non-RA")

# Build final GT Table with bold formatting
final_gt_table <- combined_table %>%
  gt::gt(
    rowname_col = "stat",    # row labels (e.g., "Median (IQR)", "Female", etc.)
    groupname_col = "var"      # group each variable together in bold
  ) %>%
  gt::tab_spanner(
    label = "Sample 1",
    columns = c("Sample1_RA", "Sample1_NonRA")
  ) %>%
  gt::tab_spanner(
    label = "Sample 2",
    columns = c("Sample2_RA", "Sample2_NonRA")
  ) %>%
  # Bold the RA/Non-RA column labels with sample sizes
  gt::cols_label(
    Sample1_RA = gt::md(sprintf("**RA (N=%d)**", n_sample1_RA)),
    Sample1_NonRA = gt::md(sprintf("**Non-RA (N=%d)**", n_sample1_NonRA)),
    Sample2_RA = gt::md(sprintf("**RA (N=%d)**", n_sample2_RA)),
    Sample2_NonRA = gt::md(sprintf("**Non-RA (N=%d)**", n_sample2_NonRA))
  ) %>%
  gt::tab_header(
    title = gt::md("**Table 1: Descriptive Summary**")
  ) %>%
  # Bold the group labels (variable names)
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups(groups = gt::everything())
  ) %>%
  # Add lines and center columns
  gt::opt_table_lines(extent = "all") %>%
  gt::cols_align("center")

# Display the final table
final_gt_table


## Table 

load("mcmc/mcmc_1b.RData")
load("mcmc/mcmc_2.RData")
load("mcmc/mcmc_3b.RData")
load("mcmc/mcmc_4.RData")

kappa.names <- c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]", "kappa[5]", "kappa[6]")

kappa_1 <- mcmc_1b[[1]][,kappa.names]
kappa_2 <- mcmc_2[[1]][,kappa.names]
kappa_3 <- mcmc_3b[[1]][,kappa.names]
kappa_4 <- mcmc_4[[1]][,kappa.names]

# Posterior Probabilities
idx <- 1:6
grid.tmp <- expand.grid(idx, idx)
grid.idx <- grid.tmp[grid.tmp[,1] < grid.tmp[,2],]
post.prob.mat <- matrix("", nrow = nrow(grid.idx), ncol = 5)
colnames(post.prob.mat) <- c("Comparison", "Model (1)", "Model (2)", "Model (3)", "Model (4)")

for (k in 1:nrow(grid.idx)) {
  
  post.prob.mat[k, 1] <- paste("Pr{", grid.idx[k,1], "<", grid.idx[k,2], "|y}", sep = "")
  post.prob.mat[k, 2] <- mean(as.numeric(kappa_1[,grid.idx[k,1]] < kappa_1[,grid.idx[k,2]]))
  post.prob.mat[k, 3] <- mean(as.numeric(kappa_2[,grid.idx[k,1]] < kappa_2[,grid.idx[k,2]]))
  post.prob.mat[k, 4] <- mean(as.numeric(kappa_3[,grid.idx[k,1]] < kappa_3[,grid.idx[k,2]]))
  post.prob.mat[k, 5] <- mean(as.numeric(kappa_4[,grid.idx[k,1]] < kappa_4[,grid.idx[k,2]]))
  
}

write.csv(post.prob.mat, file = "data/table-2.csv")

## Table 3

kappa.tbl <- matrix("", nrow = 6, ncol = 20)

kmcse.tmp1 <- apply(kappa_1, 2, mcse)
kmcse.tmp2 <- apply(kappa_2, 2, mcse)
kmcse.tmp3 <- apply(kappa_3, 2, mcse)
kmcse.tmp4 <- apply(kappa_4, 2, mcse)

for(j in 1:nrow(kappa.tbl)) {
  
  kappa.tbl[j,2] <- as.character(round(kmcse.tmp1[[j]][[1]], 3))
  kappa.tbl[j,7] <- as.character(round(kmcse.tmp2[[j]][[1]], 3))
  kappa.tbl[j,12] <- as.character(round(kmcse.tmp3[[j]][[1]], 3))
  kappa.tbl[j,17] <- as.character(round(kmcse.tmp4[[j]][[1]], 3))
  
  kappa.tbl[j,3] <- as.character(round(kmcse.tmp1[[j]][[2]], 3))
  kappa.tbl[j,8] <- as.character(round(kmcse.tmp2[[j]][[2]], 3))
  kappa.tbl[j,13] <- as.character(round(kmcse.tmp3[[j]][[2]], 3))
  kappa.tbl[j,18] <- as.character(round(kmcse.tmp4[[j]][[2]], 3))
  
}

kappa.tbl[,4] <- as.character(round(apply(kappa_1, 2, sd), 3))
kappa.tbl[,9] <- as.character(round(apply(kappa_2, 2, sd), 3))
kappa.tbl[,14] <- as.character(round(apply(kappa_3, 2, sd), 3))
kappa.tbl[,19] <- as.character(round(apply(kappa_4, 2, sd), 3))

khpd.int1 <- t(apply(kappa_1, 2, hpd))
khpd.int2 <- t(apply(kappa_2, 2, hpd))
khpd.int3 <- t(apply(kappa_3, 2, hpd))
khpd.int4 <- t(apply(kappa_4, 2, hpd))

for (i in 1:nrow(kappa.tbl)) {
  
  kappa.tbl[i,5] <- paste("(", round(khpd.int1[i,1], 3), ", ", 
                       round(khpd.int1[i,2], 3), ")", sep = "")
  
  kappa.tbl[i,10] <- paste("(", round(khpd.int2[i,1], 3), ", ", 
                        round(khpd.int2[i,2], 3), ")", sep = "")
  
  kappa.tbl[i,15] <- paste("(", round(khpd.int3[i,1], 3), ", ", 
                        round(khpd.int3[i,2], 3), ")", sep = "")
  
  kappa.tbl[i,20] <- paste("(", round(khpd.int4[i,1], 3), ", ", 
                        round(khpd.int4[i,2], 3), ")", sep = "")
  
}

kappa.tbl[,1] <- colnames(kappa_1)
kappa.tbl[,6] <- colnames(kappa_2)
kappa.tbl[,11] <- colnames(kappa_3)
kappa.tbl[,16] <- colnames(kappa_4)

colnames(kappa.tbl) <- rep(c("Parameter", "Posterior Mean", "MC Std. Err.", "Std. Dev.", "95% HPDI"), 4)

write.csv(kappa.tbl, file = "data/table-3.csv")

## Table 4

fix_eff_1 <- mcmc_1b[[1]][,-which(colnames(mcmc_1[[1]]) %in% kappa.names)]
fix_eff_2 <- mcmc_2[[1]][,-which(colnames(mcmc_2[[1]]) %in% kappa.names)]
fix_eff_3 <- mcmc_3b[[1]][,-which(colnames(mcmc_3[[1]]) %in% kappa.names)]
fix_eff_4 <- mcmc_4[[1]][,-which(colnames(mcmc_4[[1]]) %in% kappa.names)]

# Fixed Effects Table

fe.tbl <- matrix("", nrow = ncol(fix_eff_1), ncol = 20)

fix_eff_1 <- fix_eff_1[,c(31, 1:5, 32, 6:10, 33, 11:15, 34, 16:20,
                          35, 21:25, 36, 26:30)]
fix_eff_2 <- fix_eff_2[,c(19:21, 1:3, 22:24, 4:6, 25:27, 7:9, 28:30,
                          10:12, 31:33, 13:15, 34:36, 16:18)]
fix_eff_3 <- fix_eff_3[,c(31, 1:5, 32, 6:10, 33, 11:15, 34, 16:20,
                          35, 21:25, 36, 26:30)]
fix_eff_4 <- fix_eff_4[,c(19:21, 1:3, 22:24, 4:6, 25:27, 7:9, 28:30,
                          10:12, 31:33, 13:15, 34:36, 16:18)]

mcse.tmp1 <- apply(fix_eff_1, 2, mcse)
mcse.tmp2 <- apply(fix_eff_2, 2, mcse)
mcse.tmp3 <- apply(fix_eff_3, 2, mcse)
mcse.tmp4 <- apply(fix_eff_4, 2, mcse)

for(j in 1:nrow(fe.tbl)) {
  
  fe.tbl[j,2] <- as.character(round(mcse.tmp1[[j]][[1]], 3))
  fe.tbl[j,7] <- as.character(round(mcse.tmp2[[j]][[1]], 3))
  fe.tbl[j,12] <- as.character(round(mcse.tmp3[[j]][[1]], 3))
  fe.tbl[j,17] <- as.character(round(mcse.tmp4[[j]][[1]], 3))
  
  fe.tbl[j,3] <- as.character(round(mcse.tmp1[[j]][[2]], 3))
  fe.tbl[j,8] <- as.character(round(mcse.tmp2[[j]][[2]], 3))
  fe.tbl[j,13] <- as.character(round(mcse.tmp3[[j]][[2]], 3))
  fe.tbl[j,18] <- as.character(round(mcse.tmp4[[j]][[2]], 3))
  
}

fe.tbl[,4] <- as.character(round(apply(fix_eff_1, 2, sd), 3))
fe.tbl[,9] <- as.character(round(apply(fix_eff_2, 2, sd), 3))
fe.tbl[,14] <- as.character(round(apply(fix_eff_3, 2, sd), 3))
fe.tbl[,19] <- as.character(round(apply(fix_eff_4, 2, sd), 3))

hpd.int1 <- t(apply(fix_eff_1, 2, hpd))
hpd.int2 <- t(apply(fix_eff_2, 2, hpd))
hpd.int3 <- t(apply(fix_eff_3, 2, hpd))
hpd.int4 <- t(apply(fix_eff_4, 2, hpd))

for (i in 1:nrow(fe.tbl)) {
  
  fe.tbl[i,5] <- paste("(", round(hpd.int1[i,1], 3), ", ", 
                       round(hpd.int1[i,2], 3), ")", sep = "")
  
  fe.tbl[i,10] <- paste("(", round(hpd.int2[i,1], 3), ", ", 
                        round(hpd.int2[i,2], 3), ")", sep = "")
  
  fe.tbl[i,15] <- paste("(", round(hpd.int3[i,1], 3), ", ", 
                        round(hpd.int3[i,2], 3), ")", sep = "")
  
  fe.tbl[i,20] <- paste("(", round(hpd.int4[i,1], 3), ", ", 
                        round(hpd.int4[i,2], 3), ")", sep = "")
  
}

fe.tbl[,1] <- colnames(fix_eff_1)
fe.tbl[,6] <- colnames(fix_eff_2)
fe.tbl[,11] <- colnames(fix_eff_3)
fe.tbl[,16] <- colnames(fix_eff_4)

colnames(fe.tbl) <- rep(c("Parameter", "Posterior Mean", "MC Std. Err.", "Std. Dev.", "95% HPDI"), 4)

write.csv(fe.tbl, file = "data/table-4.csv")

#------------------------------------------------#
#
# change-point+ magnitude(90% prob being >0)
#
#------------------------------------------------#

onedrive<- get_business_onedrive()
file_path <- "Attachments/mcmc_trunc.RData"
temp_file <- tempfile(fileext = ".RData")
onedrive$download_file(
  src = file_path,
  dest = temp_file,
  overwrite = TRUE
)
mcmc_trunc<- get(load(temp_file)[1])
mcmc<-mcmc_trunc[[1]]

# Number of biomarkers being analyzed
k <- 6

# Define the time grid
time_grid <- seq(-20, 10, by = 0.01)

# Extract kappa and gamma from the mcmc data
kappa <- mcmc[, 7:12]
gamma <- mcmc[, 1:6]

# Define labels for each dimension
biomarker_labels <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")
time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(kappa)))

# Create a 3-dimensional array with named dimensions:
# Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
res <- array(NA, dim = c(k, length(time_grid), nrow(kappa)),
             dimnames = list(
               biomarker = biomarker_labels,
               time = time_labels,
               iteration = iteration_labels
             ))

# Loop over biomarkers, time grid, and iterations to fill the array.
for (b in 1:k) {
  gamma_b <- gamma[, b]     # gamma values for biomarker b
  kappa_b <- kappa[, b]     # kappa values for biomarker b
  
  for (i in seq_along(time_grid)) {
    t <- time_grid[i]
    for (j in 1:nrow(kappa)) {
      res[b, i, j] <- (t - kappa_b[j]) * gamma_b[j]
    }
  }
}

# Convert the array to a data frame.
result_df <- as.data.frame.table(res, responseName = "value")
colnames(result_df) <- c("biomarker", "time", "iteration", "value")

result_df$time <-as.numeric(as.character(result_df$time))

# Calculate the proportion of iterations for which 'value' > 0, for each biomarker and time point.
result_df <- result_df %>% 
  dplyr::group_by(biomarker, time) %>% 
  dplyr::summarise(prop_positive = mean(value > 0),.groups = "drop") %>% 
  arrange(desc(prop_positive))

# select time where absolute difference between prop_positive and 0.9 is minimized

closest_threshold <- result_df %>% 
  dplyr::group_by(biomarker) %>% 
  dplyr::slice(which.min(abs(prop_positive - 0.9))) %>% 
  dplyr::ungroup()%>% 
  arrange(time)



