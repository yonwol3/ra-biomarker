library(readxl)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
source("~/Github/ra-biomarker/hpd.R")

# Change-point densities + composite-time summary tables.

setwd("~/Documents/RA-Biomarker/")

model_tag <- "censcorr"                       # "censcorr" | "cens" | "trunc"
png_prefix <- c(censcorr = "censcorr_change-point-dens",
                cens     = "cens_change-point-dens",
                trunc    = "trunc_change-point-dens")[[model_tag]]
csv_prefix <- c(censcorr = "censcorr_summary",
                cens     = "censored_summary",
                trunc    = "truncated_summary")[[model_tag]]

time_grid <- seq(-20, 10, by = 0.01)

#--------------------------------------#
# Original biomarkers (Sample A)
#--------------------------------------#

load(sprintf("mcmc/mcmc_%s_A.RData", model_tag))
mcmc_A <- get(sprintf("mcmc_%s_A", model_tag))
K <- 6

delta <- mcmc_A[, 7:12]
gamma <- mcmc_A[, 1:6]

biomarkers <- c("RF IgA","RF IgM","RF IgG","ACPA IgA","ACPA IgM","ACPA IgG")
biomarker_labels <- biomarkers

# Changepoint density plots

outcome_colors <- brewer.pal(6, "Set1")
outcome_colors[6] <- "#F781BF"

png(sprintf("figures/%s_A.png", png_prefix),
    width = 1000, height = 1000, res = 100, units = "px")

plot(density(delta[,1], bw = 0.5),
     lwd = 2, col = outcome_colors[1],
     ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
     ylim = c(0, 1), xlim = c(-20, 10),
     main = "Change Point Densities (Sample A)")

for (i in 2:6) lines(density(delta[,i], bw = 0.5), lwd = 2, col = outcome_colors[i])

abline(v = 0, lty = 2, col = "blue"); abline(h = 0, lty = 1, col = "black"); grid()
legend("topleft", legend = biomarkers, col = outcome_colors, lwd = rep(2, 6), cex = 1)
dev.off()

# Table showing delta, gamma, and the composite statistic

delta_summ <- data.frame(biomarker = character(0), delta_summ = numeric(0))

for (i in 1:ncol(delta)) {
  delta_tmp <- delta[, i]
  mean <- round(mean(delta_tmp),2)
  q_2 <- round(hpd(delta_tmp)[1], 2)
  q_97 <- round(hpd(delta_tmp)[2], 2)
  delta_summ[i, 2] <- paste0(mean, "[",q_2,", ",q_97,"]")
  delta_summ[i, 1] <- biomarker_labels[i]
}

colnames(delta_summ)<- c("biomarker", "delta mean [95% HPD CrI]")

time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(delta)))

res <- array(NA, dim = c(K, length(time_grid), nrow(delta)),
             dimnames = list(biomarker = biomarker_labels, time = time_labels, iteration = iteration_labels))
for (b in 1:K) {
  gamma_tmp <- gamma[, b]; delta_tmp <- delta[, b]
  for (i in seq_along(time_grid)) {
    t <- time_grid[i]
    for (j in 1:nrow(delta)) res[b, i, j] <- max(0, t - delta_tmp[j]) * gamma_tmp[j]
  }
}

result_df <- as.data.frame.table(res, responseName = "value")
colnames(result_df) <- c("biomarker", "time", "iteration", "value")
result_df$time <- as.numeric(as.character(result_df$time))

result_df <- result_df %>%
  dplyr::group_by(biomarker, time) %>%
  dplyr::summarise(prop_positive = mean(value > 0), .groups = "drop") %>%
  arrange(desc(prop_positive))

closest_threshold <- result_df %>%
  dplyr::group_by(biomarker) %>%
  dplyr::slice(which.min(abs(prop_positive - 0.9))) %>%
  dplyr::ungroup() %>% arrange(time)

gamma_summ <- matrix(NA, ncol = 2, nrow = length(biomarker_labels))

for (i in 1:length(biomarker_labels)) {
  gamma_tmp <- gamma[, i]
  mean <- round(mean(gamma_tmp),2)
  q_2 <- round(hpd(gamma_tmp)[1], 2)
  q_97 <- round(hpd(gamma_tmp)[2], 2)
  gamma_summ[i, 2] <- paste0(mean, "[",q_2,", ",q_97,"]")
  gamma_summ[i, 1] <- biomarker_labels[i]
}

gamma_summ <- as.data.frame(gamma_summ)
colnames(gamma_summ) <- c("biomarker", "gamma mean [95% HPD CrI]")
closest_threshold <- left_join(closest_threshold, delta_summ, by = "biomarker")
closest_threshold <- left_join(closest_threshold, gamma_summ, by = "biomarker")
write.csv(closest_threshold, sprintf("tables/%s_A.csv", csv_prefix))

#--------------------------------------#
# New biomarkers (Sample B)
#--------------------------------------#

load(sprintf("mcmc/mcmc_%s_B.RData", model_tag))
mcmc_B <- get(sprintf("mcmc_%s_B", model_tag))
K <- 8

delta <- mcmc_B[, 9:16]
gamma <- mcmc_B[, 1:8]

time_labels <- as.character(time_grid)
iteration_labels <- paste0("iter", seq_len(nrow(delta)))

biomarkers <- c("anti-CCP3 (IgG)","anti-citVim2 (IgG)", "anti-citFib (IgG)","anti-citHis1 (IgG)",
                "anti-CCP3 (IgA)","anti-citVim2 (IgA)","anti-citFib (IgA)","anti-citHis1 (IgA)")
biomarker_labels <- biomarkers

outcome_colors <- brewer.pal(8, "Paired")
names(outcome_colors) <- biomarkers

# Changepoint density plots

png(sprintf("figures/%s_B.png", png_prefix),
    width = 1000, height = 1000, res = 100, units = "px")

plot(density(delta[,1], bw = 0.5),
     lwd = 2, col = outcome_colors[1],
     ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
     ylim = c(0, 0.8), xlim = c(-20, 10),
     main = "Change Point Densities (Sample B)")
for (i in 2:8) lines(density(delta[,i], bw = 0.5), lwd = 2, col = outcome_colors[i])
abline(v = 0, lty = 2, col = "blue"); abline(h = 0, lty = 1, col = "black"); grid()
legend("topleft", legend = biomarker_labels, col = outcome_colors, lwd = rep(2, 8), cex = 1)
dev.off()

# Table showing delta, gamma, and the composite statistic

delta_summ <- data.frame(biomarker = character(0), delta_summ = numeric(0))

for (i in 1:ncol(delta)) {
  delta_tmp <- delta[, i]
  mean <- round(mean(delta_tmp),2)
  q_2 <- round(hpd(delta_tmp)[1], 2)
  q_97 <- round(hpd(delta_tmp)[2], 2)
  delta_summ[i, 2] <- paste0(mean, "[",q_2,", ",q_97,"]")
  delta_summ[i, 1] <- biomarker_labels[i]
}

colnames(delta_summ) <- c("biomarker", "delta mean [95% HPD CrI]")

res <- array(NA, dim = c(K, length(time_grid), nrow(delta)),
             dimnames = list(biomarker = biomarker_labels, time = time_labels, iteration = iteration_labels))

for (b in 1:K) {
  gamma_tmp <- gamma[, b]; delta_tmp <- delta[, b]
  for (i in seq_along(time_grid)) {
    t <- time_grid[i]
    for (j in 1:nrow(delta)) res[b, i, j] <- max(0, t - delta_tmp[j]) * gamma_tmp[j]
  }
}

result_df <- as.data.frame.table(res, responseName = "value")
colnames(result_df) <- c("biomarker", "time", "iteration", "value")
result_df$time <- as.numeric(as.character(result_df$time))

result_df <- result_df %>%
  dplyr::group_by(biomarker, time) %>%
  dplyr::summarise(prop_positive = mean(value > 0), .groups = "drop") %>%
  arrange(desc(prop_positive))

closest_threshold <- result_df %>%
  dplyr::group_by(biomarker) %>%
  dplyr::slice(which.min(abs(prop_positive - 0.9))) %>%
  dplyr::ungroup() %>% arrange(time)

gamma_summ <- matrix(NA, ncol = 2, nrow = length(biomarker_labels))

for (i in 1:length(biomarker_labels)) {
  gamma_tmp <- gamma[, i]
  mean <- round(mean(gamma_tmp),2)
  q_2 <- round(hpd(gamma_tmp)[1], 2)
  q_97 <- round(hpd(gamma_tmp)[2], 2)
  gamma_summ[i, 2] <- paste0(mean, "[",q_2,", ",q_97,"]")
  gamma_summ[i, 1] <- biomarker_labels[i]
}

gamma_summ <- as.data.frame(gamma_summ)
colnames(gamma_summ)<- c("biomarker","gamma mean [95% HPD CrI]")
closest_threshold <- left_join(closest_threshold, delta_summ, by = "biomarker")
closest_threshold <- left_join(closest_threshold, gamma_summ, by = "biomarker")
write.csv(closest_threshold, sprintf("tables/%s_B.csv", csv_prefix))

#--------------------------------------#
# Main-text Figure 3: combined (A | B)
#--------------------------------------#

if (model_tag == "censcorr") {

  library(magick)
  library(cowplot)

  img1 <- image_read(sprintf("figures/%s_A.png", png_prefix))
  img2 <- image_read(sprintf("figures/%s_B.png", png_prefix))

  g1 <- ggdraw() +
    draw_image(img1, scale = 1) +
    theme(plot.margin = unit(rep(0, 4), "cm"))

  g2 <- ggdraw() +
    draw_image(img2, scale = 1) +
    theme(plot.margin = unit(rep(0, 4), "cm"))

  grid_plot <- plot_grid(
    g1, g2,
    labels = c("A", "B"),
    label_size = 22,
    ncol = 2,
    align = "hv"
  )

  png("figures/fig3.png", height = 800, width = 1000)
  print(grid_plot)
  dev.off()

}
