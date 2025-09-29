library(readxl)
library(tidyverse)
library(Microsoft365R)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(kableExtra)
library(magick)
library(cowplot)
source("~/Github/ra-biomarker/hpd.R")

#--------------------------------------#
# Original biomarkers w/truncation at LOD
#--------------------------------------#

setwd("~/Documents/RA-Biomarker/")
temp_file <- "mcmc/mcmc_trunc_A.RData"
time_grid <- seq(-20, 10, by = 0.01)
K <- 6

mcmc <- load(temp_file)
kappa <- mcmc_trunc_A[ ,7:12]
gamma <- mcmc_trunc_A[ ,1:6]

biomarkers <- c("RF IgA","RF IgM","RF IgG","ACPA IgA","ACPA IgM","ACPA IgG")
biomarker_labels <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")

### changepoint density plots####
outcome_colors <- brewer.pal(6, "Set1")
outcome_colors[6] <- "#F781BF"

png("figures/trunc_change-point-dens_A.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(kappa[,1]), lwd = 2,
     col = outcome_colors[1], ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
     ylim = c(0, 1.5),
     xlim = c(-20, 5),
     main = "Change Point Densities (Sample A)")

for (i in 2:6) {
  lines(density(kappa[,i]), lwd = 2, col = outcome_colors[i])
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

### Table showing kappa, gamma, and the statistic (gamma + kappa)

kappa_summ <- data.frame(biomarker = character(0), kappa_summ = numeric(0)) 

for (i in 1:ncol(kappa)) {
  kappa_tmp <- kappa[, i]
  mean <- round(mean(kappa_tmp),2)
  q_2 <- round(hpd(kappa_tmp)[1], 2) # 2.5th quantile
  q_97 <-round(hpd(kappa_tmp)[2], 2) 
  kappa_summ[i, 2] <- paste(mean, "[",q_2,"-",q_97,"]")
  kappa_summ[i, 1] <- biomarker_labels[i]
}

colnames(kappa_summ) <- c("biomarker", "kappa mean[95% HPD CrI]")

# Define labels for each dimension

time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(kappa)))

# Create a 3-dimensional array with named dimensions:
# Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
res <- array(NA, dim = c(K, length(time_grid), nrow(kappa)),
             dimnames = list(
               biomarker = biomarker_labels,
               time = time_labels,
               iteration = iteration_labels)
             )

# Loop over biomarkers, time grid, and iterations to fill the array.
for (b in 1:K) {
  
  gamma_tmp <- gamma[, b]     # gamma values for biomarker b
  kappa_tmp <- kappa[, b]     # kappa values for biomarker b
  
  for (i in seq_along(time_grid)) {
    
    t <- time_grid[i]
    
    for (j in 1:nrow(kappa)) {
      res[b, i, j] <- (t - kappa_tmp[j]) * gamma_tmp[j]
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

gamma_summ <- matrix(NA,ncol = 2, nrow=length(biomarker_labels))

for (i in 1:length(biomarker_labels)) {
  gamma_b <- gamma[, i]
  mean <- round(mean(gamma_b),2)
  q_2 <- round(hpd(gamma_b)[1], 2) # 2.5th quantile
  q_97 <- round(hpd(gamma_b)[2], 2) # 97.5th quantile
  gamma_summ[i, 2] <- paste(mean, "[",q_2,"-",q_97,"]")
  gamma_summ[i, 1] <- biomarker_labels[i]
}

gamma_summ <- as.data.frame(gamma_summ)
colnames(gamma_summ) <- c("biomarker","gamma mean[95% HPD CrI]")
closest_threshold <- left_join(closest_threshold,kappa_summ, by = "biomarker") 
closest_threshold <- left_join(closest_threshold, gamma_summ, by = "biomarker")
write.csv(closest_threshold,"tables/truncated_summary_A.csv")

#--------------------------------------#
# New biomarkers w/truncation at LOD
#--------------------------------------#

load("mcmc/mcmc_trunc_B.RData")

kappa <- mcmc_trunc_B[, 9:16]
gamma <- mcmc_trunc_B[, 1:8]

# Define labels for each dimension
biomarkers<-c("anti-CCP3 (IgG)","anti-citVim2 (IgG)", "anti-citFib (IgG)","anti-citHis1 (IgG)",
              "anti-CCP3 (IgA)","anti-citVim2 (IgA)","anti-citFib (IgA)","anti-citHis1 (IgA)")
biomarker_labels <- biomarkers

outcome_colors <- brewer.pal(8, "Paired")
names(outcome_colors) <- biomarkers

### changepoint density plots####

png("figures/trunc_change-point-dens_B.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(kappa[,1]), lwd = 2,
     col = outcome_colors[1], ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
     ylim = c(0, 0.8),
     xlim = c(-20, 5),
     main = "Change Point Densities (Sample B)")

for (i in 2:8) {
  lines(density(kappa[,i]), lwd = 2, col = outcome_colors[i])
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

### Table showing kappa, gamma, and the statistic (gamma + kappa)

kappa_summ <- data.frame(biomarker = character(0),kappa_summ = numeric(0))

for (i in 1:ncol(kappa)) {
  kappa_tmp <- kappa[, i]
  mean <- round(mean(kappa_tmp),2)
  q_2 <- round(hpd(kappa_tmp)[1], 2) # 2.5th quantile
  q_97 <- round(hpd(kappa_tmp)[2], 2) 
  kappa_summ[i, 2] <- paste(mean, "[",q_2,"-",q_97,"]")
  kappa_summ[i, 1] <- biomarker_labels[i]
}

colnames(kappa_summ) <- c("biomarker","kappa mean[95% HPD CrI]")

K <- 8
time_grid <- seq(-20, 10, by = 0.01)
time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(kappa)))

# Create a 3-dimensional array with named dimensions:
# Dimension 1: Biomarker, Dimension 2: Time, Dimension 3: Iteration
res <- array(NA, dim = c(K, length(time_grid), nrow(kappa)),
             dimnames = list(
               biomarker = biomarker_labels,
               time = time_labels,
               iteration = iteration_labels)
             )

# Loop over biomarkers, time grid, and iterations to fill the array.
for (b in 1:K) {
  
  gamma_tmp <- gamma[, b]     # gamma values for biomarker b
  kappa_tmp <- kappa[, b]     # kappa values for biomarker b
  
  for (i in seq_along(time_grid)) {
    
    t <- time_grid[i]
    
    for (j in 1:nrow(kappa)) {
      res[b, i, j] <- (t - kappa_tmp[j]) * gamma_tmp[j]
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

for ( i in 1:length(biomarker_labels)) {
  gamma_tmp <- gamma[, i]
  mean <- round(mean(gamma_tmp),2)
  q_2 <- round(hpd(gamma_tmp)[1], 2) # 2.5th quantile
  q_97 <- round(hpd(gamma_tmp)[2], 2) # 97.5th quantile
  gamma_summ[i, 2] <- paste(mean, "[",q_2,"-",q_97,"]")
  gamma_summ[i, 1] <- biomarker_labels[i]
}

gamma_summ <- as.data.frame(gamma_summ)
colnames(gamma_summ) <- c("biomarker","gamma mean[95% HPD CrI]")
closest_threshold <- left_join(closest_threshold,kappa_summ, by="biomarker") 
closest_threshold <- left_join(closest_threshold, gamma_summ, by="biomarker")
write.csv(closest_threshold,"tables/truncated_summary_B.csv")

# combining the density plots into one plot for the main figure
img1 <- image_read("figures/trunc_change-point-dens_A.png")
img2 <- image_read("figures/trunc_change-point-dens_B.png")

# Create ggplot-like objects from the images
g1 <- ggdraw() +
  draw_image(img1, scale = 1) +
  theme(plot.margin = unit(rep(0, 4), "cm"))

g2 <- ggdraw() +
  draw_image(img2, scale = 1) +
  theme(plot.margin = unit(rep(0, 4), "cm"))

# Arrange side by side with labels A and B
grid_plot <- plot_grid(
  g1, g2,
  labels     = c("A", "B"),
  label_size = 22,
  ncol       = 2,
  align      = "hv"
)

png("~/Documents/RA-Biomarker/figures/fig3.png", height = 800, width = 1000)
print(grid_plot)
dev.off()
