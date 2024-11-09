##########################################################
## PURPOSE: Script for creation of scatterplots and     ##
##          smoothed splines for exploratory analysis   ##
##          of RA dataset                               ##
## BY:      Kevin Josey                                 ##
##########################################################

# Libraries

library(tidyverse)
library(rjags)
library(plyr)
library(splines)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# Reading data cleaning R file
source("clean-data.R")



# Naming serum antibodies correctly
colnames(Y) <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")
clean<- clean %>%
         dplyr::rename('RF IgA'= igarfconc_ , 
                'RF IgM'= igmrfconc_ , 
                'RF IgG'= iggrfconc_ , 
                'ACPA IgA'= igaccpavgconc, 
                'ACPA IgM'= igmccpavgconc, 
                'ACPA IgG'= iggccpavgconc)


## Descriptive Spline Plot

png(filename = "figures/spline-plots-horiz.png",
    width = 12,
    height = 8,
    units = "in",
    res = 300)

par(mfrow = c(2,3))

# smoothing splines - igarfconc_
sm_spline1 <- smooth.spline(raDat$t_yrs, Y[,1], df = 4)
plot(raDat$t_yrs, Y[,1], col="grey",
     xlab = "Time before Diagnosis", ylab = "RF-IgA",
     main = "RF-IgA Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline1, col = "darkolivegreen4", lwd = 2)

#smoothing splines - igmrfconc_
sm_spline2 <- smooth.spline(raDat$t_yrs, Y[,2], df = 4)
plot(raDat$t_yrs, Y[,2], col = "grey",
     xlab = "Time before Diagnosis", ylab = "RF-IgM", 
     main = "RF-IgM Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline2, col = "darkorange4", lwd = 2)

#smoothing splines - iggrfconc_
sm_spline3 <- smooth.spline(raDat$t_yrs, Y[,3], df = 4)
plot(raDat$t_yrs, Y[,3], col = "grey",
     xlab = "Time before Diagnosis", ylab="RF-IgG", 
     main = "RF-IgG Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline3, col= "firebrick4", lwd = 2)

#smoothing splines - igaccpavgconc
sm_spline4 <- smooth.spline(raDat$t_yrs, Y[,4], df = 4)
plot(raDat$t_yrs, Y[,4],col = "grey",
     xlab = "Time before Diagnosis", ylab="ACPA-IgA",
     main="ACPA-IgA Serum Levels over Time with \nSmoothing Spline")
abline(v=0, lty=3, col = "black")
lines(sm_spline4, col = "dodgerblue4", lwd = 2)

#smoothing splines - igmccpavgconc
sm_spline5 <- smooth.spline(raDat$t_yrs, Y[,5], df = 4)
plot(raDat$t_yrs, Y[,5], col="grey",
     xlab = "Time before Diagnosis", ylab = "ACPA-IgM",
     main = "ACPA-IgM Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline5, col = "mediumorchid4", lwd = 2)

sm_spline6 <- smooth.spline(raDat$t_yrs, Y[,6], df = 4)
plot(raDat$t_yrs, Y[,6], col="grey",
     xlab = "Time before Diagnosis", ylab = "ACPA-IgG",
     main = "ACPA-IgG Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline6, col= "khaki4", lwd = 2)
dev.off()


## Changepoint Density Plot

load("mcmc/mcmc_c.RData")

kappa.names <- c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]", "kappa[5]", "kappa[6]")

kappa <- mcmc_c[[1]][,kappa.names]
colnames(kappa) <- colnames(Y)

png("figures/change-point-dist.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(kappa[,1]), lwd = 2,
     col = "darkolivegreen4", ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
     ylim = c(0, 0.8),
     xlim = c(-20, 5),
     main = "Change Point Densities")
lines(density(kappa[,2]), lwd = 2,
      col = "darkorange4")
lines(density(kappa[,3]), lwd = 2,
      col = "firebrick4")
lines(density(kappa[,4]), lwd = 2,
      col = "dodgerblue4")
lines(density(kappa[,5]), lwd = 2,
      col = "mediumorchid4")
lines(density(kappa[,6]), lwd = 2,
      col = "khaki4")
abline(v = 0, lty = 2, col = "blue")
abline(h = 0, lty = 1, col = "black")
grid()
legend("topleft", 
       legend = colnames(Y),
       col = c("darkolivegreen4", "darkorange4", "firebrick4", 
               "dodgerblue4", "mediumorchid4", "khaki4"), 
       lwd = c(2, 2, 2, 2, 2, 2),
       cex = 1)

dev.off()

## Diagnostics of model (c)

load("mcmc/mcmc_c.RData")

pdf("images/trace_c.pdf")
plot(mcmc_c)
dev.off()

mcmc <- as.matrix(mcmc_c[[1]])

png("images/acf_c.png", 
    width = 7000, 
    height = 6000,
    res = 100, 
    units = "px")
par(mfrow = c(4, 6), mar = c(1,1,1,1))
for (i in 1:ncol(mcmc))
  acf(mcmc[,i], ylab = colnames(mcmc)[i])
title("ACF Plots", outer = TRUE)
dev.off()

## Diagnostics of model (c)

load("mcmc/mcmc_c.RData")

pdf("images/trace_c.pdf")
plot(mcmc_c)
dev.off()

mcmc <- as.matrix(mcmc_c[[1]])

png("images/acf_c.png", 
    width = 7000, 
    height = 6000,
    res = 100, 
    units = "px")
par(mfrow = c(4, 3), mar = c(1,1,1,1))
for (i in 1:ncol(mcmc))
  acf(mcmc[,i], ylab = colnames(mcmc)[i])
title("ACF Plots", outer = TRUE)
dev.off()






#################################################
## Visualize biomarkers for different RA groups
#################################################

# Assuming the original 'diagnosis' variable uses 0 for 'no_RA' and 1 for 'RA'
clean$diagnosis <- factor(clean$diagnosis, levels = c(0, 1), labels = c("no_RA", "RA"))

# List of outcome variables
outcome_vars <-colnames(Y)

# Color palette for outcomes (Set1 is suitable for up to 9 distinct colors)
outcome_colors <- brewer.pal(6, "Set1")
names(outcome_colors) <- outcome_vars

# Convert 'diagnosis' to a factor with explicit levels and labels

# Define colors for 'RA' and 'no_RA'
group_colors <- c("RA" = "dodgerblue4", "no_RA" = "grey")  # Blue for RA, gray for no_RA

# Transparency for scatter points
point_alpha <- 0.5

# Initialize a list to store individual plots
plot_list <- list()

# Loop through each outcome variable to create plots
for (i in seq_along(outcome_vars)) {
  outcome_var <- outcome_vars[i]
  color <- outcome_colors[outcome_var]
  
  # Prepare data for plotting
  df_plot <- clean[, c("time", "diagnosis", outcome_var)]
  names(df_plot)[3] <- "outcome"  # Rename the outcome column for consistency
  
  # Fit smooth.spline for 'RA'
  df_ra <- df_plot[df_plot$diagnosis == "RA", ]
  spline_ra <- smooth.spline(df_ra$time, df_ra$outcome, df = 4)
  
  # Fit smooth.spline for 'no_RA'
  df_no_ra <- df_plot[df_plot$diagnosis == "no_RA", ]
  spline_no_ra <- smooth.spline(df_no_ra$time, df_no_ra$outcome, df = 4)
  
  # Create a sequence of time points for smooth curves
  time_seq <- seq(min(clean$time, na.rm = TRUE), max(clean$time, na.rm = TRUE), length.out = 1000)
  
  # Predict values using the fitted splines
  pred_ra <- predict(spline_ra, x = time_seq)
  pred_no_ra <- predict(spline_no_ra, x = time_seq)
  
  # Create data frames for the fitted spline lines
  df_spline_ra <- data.frame(time = pred_ra$x, outcome = pred_ra$y, diagnosis = "RA")
  df_spline_no_ra <- data.frame(time = pred_no_ra$x, outcome = pred_no_ra$y, diagnosis = "no_RA")
  
  # Generate the plot
  p <- ggplot() +
    # Scatter points for actual data
    geom_point(data = df_plot, 
               aes(x = time, y = outcome, color = diagnosis), 
               alpha = point_alpha) +
    # Manual color scale
    scale_color_manual(values = group_colors) +
    # Smoothing spline for 'RA' (solid line)
    geom_line(data = df_spline_ra, 
              aes(x = time, y = outcome), 
              color = color, 
              size = 1) +
    # Smoothing spline for 'no_RA' (dashed line)
    geom_line(data = df_spline_no_ra, 
              aes(x = time, y = outcome), 
              color = color, 
              linetype = "dashed", 
              size = 1) +
    # Labels and title
    labs(title = paste(outcome_var, "Serum Levels over Time with Smoothing Spline"),
         x = "Time",
         y = outcome_var) +
    # Minimal theme for a clean look
    theme_minimal() +
    # Center the plot title and remove legend title
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # Add the plot to the list
  plot_list[[i]] <- p
}

  grid.arrange(grobs = plot_list, ncol = 2)





