library(readxl)
library(tidyverse)
library(Microsoft365R)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(kableExtra)
source("~/Github/ra-biomarker/hpd.R")

#--------------------------------------#
# Original biomarkers w/censoring at LOD
#--------------------------------------#

setwd("~/Documents/RA-Biomarker/")
load("mcmc/mcmc_cens_A.RData")
K <- 6
time_grid <- seq(-20, 10, by = 0.01)

phi <- mcmc_cens_A[ ,13:18] 
delta <- mcmc_cens_A[ ,7:12] 
gamma <- mcmc_cens_A[ ,1:6]

biomarkers <- c("RF IgA","RF IgM","RF IgG","ACPA IgA","ACPA IgM","ACPA IgG")
biomarker_labels <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")

### changepoint density plots####
outcome_colors <- brewer.pal(6, "Set1")
outcome_colors[6] <- "#F781BF"

png("figures/cens_change-point-dens_A.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(delta[,1]), 
     lwd = 2,
     col = outcome_colors[1], 
     ylab = "Posterior Density", 
     xlab = "Years Prior to Diagnosis",
     ylim = c(0, 1),
     xlim = c(-20, 5),
     main = "Change Point Densities (Sample A)")

for (i in 2:6) {
  lines(density(delta[,i]), lwd = 2, col = outcome_colors[i])
}

abline(v = 0, lty = 2, col = "blue")
abline(h = 0, lty = 1, col = "black")
grid()

legend("topleft", 
       legend = biomarkers,
       col = outcome_colors, 
       lwd = rep(2, 6),
       cex = 1)

dev.off()

### Table showing delta, gamma, and the statistic (gamma + delta)

delta_summ <- data.frame(biomarker = character(0), delta_summ = numeric(0))

for (i in 1:ncol(delta)) {
  delta_tmp <- delta[, i]
  mean <- round(mean(delta_tmp),2)
  q_2 <- round(hpd(delta_tmp)[1], 2) # 2.5th quantile
  q_97 <- round(hpd(delta_tmp)[2], 2) 
  delta_summ[i, 2] <- paste(mean, "[",q_2,", ",q_97,"]")
  delta_summ[i, 1] <- biomarker_labels[i]
}

colnames(delta_summ)<- c("biomarker", "delta mean [95% HPD CrI]")

# Define labels for each dimension

time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(delta)))

# Create a 3-dimensional array with named dimensions:
# Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
res <- array(NA, dim = c(K, length(time_grid), nrow(delta)),
             dimnames = list(
               biomarker = biomarker_labels,
               time = time_labels,
               iteration = iteration_labels)
             )

# Loop over biomarkers, time grid, and iterations to fill the array.
for (b in 1:K) {
  
  gamma_tmp <- gamma[, b]     # gamma values for biomarker b
  delta_tmp <- delta[, b]     # delta values for biomarker b
  phi_tmp <- phi[, b]     # phi values for biomarker b
  
  for (i in seq_along(time_grid)) {
    
    t <- time_grid[i]
    
    for (j in 1:nrow(delta)) {
      res[b, i, j] <- (t - delta_tmp[j]) * gamma_tmp[j] * exp((t - delta_tmp[j])*phi_tmp[j])
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

gamma_summ <- matrix(NA, ncol = 2, nrow = length(biomarker_labels))

for (i in 1:length(biomarker_labels)) {
  gamma_tmp <- gamma[, i]
  mean <- round(mean(gamma_tmp),2)
  q_2 <- round(hpd(gamma_tmp)[1], 2) # 2.5th quantile
  q_97 <- round(hpd(gamma_tmp)[2], 2) # 97.5th quantile
  gamma_summ[i, 2] <- paste(mean, "[",q_2,", ",q_97,"]")
  gamma_summ[i, 1] <- biomarker_labels[i]
}

gamma_summ <- as.data.frame(gamma_summ)
colnames(gamma_summ) <- c("biomarker", "gamma mean [95% HPD CrI]")
closest_threshold <- left_join(closest_threshold, delta_summ, by = "biomarker") 
closest_threshold <- left_join(closest_threshold, gamma_summ, by = "biomarker")
write.csv(closest_threshold,"tables/censored_summary_A.csv")

#--------------------------------------#
# New biomarkers w/censoring at LOD
#--------------------------------------#

load("mcmc/mcmc_cens_B.RData")
K <- 8
time_grid <- seq(-20, 10, by = 0.01)
time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(delta)))

phi <- mcmc_cens_B[, 17:24]
delta <- mcmc_cens_B[, 9:16]
gamma <- mcmc_cens_B[, 1:8]

# Define labels for each dimension
biomarkers <- c("anti-CCP3 (IgG)","anti-citVim2 (IgG)", "anti-citFib (IgG)","anti-citHis1 (IgG)",
                "anti-CCP3 (IgA)","anti-citVim2 (IgA)","anti-citFib (IgA)","anti-citHis1 (IgA)")
biomarker_labels <- biomarkers

outcome_colors <- brewer.pal(8, "Paired")
names(outcome_colors) <- biomarkers

### changepoint density plots####

png("figures/cens_change-point-dens_B.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(delta[,1]), 
     lwd = 2,
     col = outcome_colors[1], 
     ylab = "Posterior Density", 
     xlab = "Years Prior to Diagnosis",
     ylim = c(0, 0.8),
     xlim = c(-20, 5),
     main = "Change Point Densities (Sample B)")

for (i in 2:8) {
  lines(density(delta[,i]), lwd = 2, col = outcome_colors[i])
}

abline(v = 0, lty = 2, col = "blue")
abline(h = 0, lty = 1, col = "black")
grid()

legend("topleft", 
       legend = biomarker_labels,
       col = outcome_colors, 
       lwd = rep(2, 8),
       cex = 1)

dev.off()

### Table showing delta, gamma, and the statistic (gamma + delta)

delta_summ <- data.frame(biomarker = character(0), delta_summ = numeric(0))

for (i in 1:ncol(delta)) {
  delta_tmp <- delta[, i]
  mean <- round(mean(delta_tmp),2)
  q_2 <- round(hpd(delta_tmp)[1], 2) # 2.5th quantile
  q_97 <- round(hpd(delta_tmp)[2], 2) 
  delta_summ[i, 2] <- paste(mean, "[",q_2,", ",q_97,"]")
  delta_summ[i, 1] <- biomarker_labels[i]
}

colnames(delta_summ) <- c("biomarker", "delta mean [95% HPD CrI]")

# Create a 3-dimensional array with named dimensions:
# Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
res <- array(NA, dim = c(K, length(time_grid), nrow(delta)),
             dimnames = list(
               biomarker = biomarker_labels,
               time = time_labels,
               iteration = iteration_labels))

# Loop over biomarkers, time grid, and iterations to fill the array.
for (b in 1:K) {
  
  gamma_tmp <- gamma[, b]     # gamma values for biomarker b
  delta_tmp <- delta[, b]     # delta values for biomarker b
  phi_tmp <- phi[, b]     # phi values for biomarker b
  
  for (i in seq_along(time_grid)) {
    
    t <- time_grid[i]
    
    for (j in 1:nrow(delta)) {
      res[b, i, j] <- (t - delta_tmp[j]) * gamma_tmp[j] * exp((t - delta_tmp[j])*phi_tmp[j])
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

gamma_summ <- matrix(NA, ncol = 2, nrow = length(biomarker_labels))

for (i in 1:length(biomarker_labels)) {
  gamma_tmp <- gamma[, i]
  mean <- round(mean(gamma_tmp),2)
  q_2 <- round(hpd(gamma_tmp)[1], 2) # 2.5th quantile
  q_97 <-round(hpd(gamma_tmp)[2], 2) # 97.5th quantile
  gamma_summ[i, 2] <- paste(mean, "[",q_2,", ",q_97,"]")
  gamma_summ[i, 1] <- biomarker_labels[i]
}

gamma_summ <- as.data.frame(gamma_summ)
colnames(gamma_summ)<- c("biomarker","gamma mean [95% HPD CrI]")
closest_threshold<-left_join(closest_threshold, delta_summ, by = "biomarker") 
closest_threshold<-left_join(closest_threshold, gamma_summ, by = "biomarker")
write.csv(closest_threshold,"tables/censored_summary_B.csv")
