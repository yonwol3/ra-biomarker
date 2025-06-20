##########################################################
## PURPOSE: Script for creation of scatterplots and     ##
##          smoothed splines for exploratory analysis   ##
##          of RA dataset                               ##
## BY:      Yonatan Wolde                               ##
##########################################################

# Libraries

library(tidyverse)
library(plyr)
library(splines)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(kableExtra)

# Reading data cleaning R file
# setwd("~/Dropbox/Projects/RA-Biomarker/")
# source("~/Github/ra-biomarker/clean-data.R")
source("clean-data.R")

# Naming serum antibodies correctly
colnames(Y) <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")
clean <- clean %>%
  dplyr::rename('RF IgA'= igarfconc_ , 
                'RF IgM'= igmrfconc_ , 
                'RF IgG'= iggrfconc_ , 
                'ACPA IgA'= igaccpavgconc, 
                'ACPA IgM'= igmccpavgconc, 
                'ACPA IgG'= iggccpavgconc)

## LOESS Plot

# Assuming the original 'diagnosis' variable uses 0 for 'no_RA' and 1 for 'RA'
clean$diagnosis <- factor(clean$diagnosis, levels = c(0, 1), labels = c("No RA", "RA"))

# List of outcome variables
outcome_vars <-colnames(Y)

# Color palette for outcomes (Set1 is suitable for up to 9 distinct colors)
outcome_colors <- brewer.pal(6, "Set1")
outcome_colors[6] <- "#F781BF" # replace the yellow color
names(outcome_colors) <- outcome_vars

# Convert 'diagnosis' to a factor with explicit levels and labels

# Define colors for 'RA' and 'no_RA'
group_colors <- c("RA" = "dodgerblue4", "No RA" = "grey")  # Blue for RA, gray for no_RA

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
  df_no_ra <- df_plot[df_plot$diagnosis == "No RA", ]
  spline_no_ra <- smooth.spline(df_no_ra$time, df_no_ra$outcome, df = 4)
  
  # Create a sequence of time points for smooth curves
  time_seq <- seq(min(clean$time, na.rm = TRUE), max(clean$time, na.rm = TRUE), length.out = 1000)
  
  # Predict values using the fitted splines
  pred_ra <- predict(spline_ra, x = time_seq)
  pred_no_ra <- predict(spline_no_ra, x = time_seq)
  
  # Create data frames for the fitted spline lines
  df_spline_ra <- data.frame(time = pred_ra$x, outcome = pred_ra$y, diagnosis = "RA")
  df_spline_no_ra <- data.frame(time = pred_no_ra$x, outcome = pred_no_ra$y, diagnosis = "No RA")
  
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
    geom_vline(xintercept = 0, linetype="dashed") +
    # Labels and title
    labs(title = paste(outcome_var, "Serum Levels over Time with Smoothing Spline"),
         x = "Time Prior to Diagnosis ",
         y = paste(outcome_var)) +
    # Minimal theme for a clean look
    theme_minimal() +
    # Center the plot title and remove legend title
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # Add the plot to the list
  plot_list[[i]] <- p
  
}


ggarrange(plot_list[[1]],plot_list[[2]], plot_list[[3]],
          plot_list[[4]],plot_list[[5]], plot_list[[6]], nrow = 3, ncol = 2, 
          legend="right",align = "v", common.legend = TRUE)


#---------------------------------------------------#
# change-point density for the original biomarkers
#---------------------------------------------------#
# source("clean-data.R")
# biomarkers<-c("RF IgA","RF IgM","RF IgG","ACPA IgA","ACPA IgM","ACPA IgG")
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
# kappa<-mcmc[,7:12] # change points only (no gamma (rate of change))
# 
# png("figures/change-point-dens-originalbiomarkers.png", 
#     width = 1000, 
#     height = 1000,
#     res = 100, 
#     units = "px")
# 
# plot(density(kappa[,1]), lwd = 2,
#      col = outcome_colors[1], ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
#      ylim = c(0, 1.5),
#      xlim = c(-20, 5),
#      main = "Change Point Densities")
# 
# for (i in 2:6) {
#   lines(density(kappa[,i]), lwd = 2, col = outcome_colors[i])
# }
# 
# abline(v = 0, lty = 2, col = "blue")
# abline(h = 0, lty = 1, col = "black")
# grid()
# 
# legend("topleft", 
#        legend = biomarkers,
#        col = outcome_colors, 
#        lwd = rep(2, 6),
#        cex = 1)
# 
# dev.off()
# 
# #-------------------------------------------#
# # calculate mean and 95% Credible Interval 
# #-------------------------------------------#
# credible_res <- data.frame(biomarker=character(0),mean = numeric(0), 
#                            credible_lower = numeric(0),
#                            credible_upper = numeric(0))
# 
# for (i in 1:ncol(kappa)) {
#   credible_res[i, "mean"] <- round(mean(kappa[, i]), 2)         # mean 
#   credible_res[i, "credible_lower"] <- round(quantile(kappa[, i], 0.025), 2)  # 2.5th quantile
#   credible_res[i, "credible_upper"] <- round(quantile(kappa[, i], 0.975), 2)  # 97.5th quantile
#   credible_res[i, "biomarker"]<- biomarkers[i]
# }
# credible_res<-arrange(credible_res,mean)
# credible_res<-credible_res %>% 
#               mutate(`95% Credible Interval`=paste("[", credible_lower, ", ", credible_upper, "]"))%>%
#               select(biomarker, mean, `95% Credible Interval`)
#              
# write.csv(credible_res,"../../original_cp_ci.csv")


