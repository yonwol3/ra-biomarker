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

### Sample A

# Reading data cleaning R file
setwd("~/Dropbox/Projects/RA-Biomarker/")
source("~/Github/ra-biomarker/clean-data-A.R")

# Naming serum antibodies correctly
colnames(Y) <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")
data_A <- data_A %>%
  dplyr::rename('RF IgA'= igarfconc_ , 
                'RF IgM'= igmrfconc_ , 
                'RF IgG'= iggrfconc_ , 
                'ACPA IgA'= igaccpavgconc, 
                'ACPA IgM'= igmccpavgconc, 
                'ACPA IgG'= iggccpavgconc)

# Assuming the original 'diagnosis' variable uses 0 for 'no_RA' and 1 for 'RA'
data_A$diagnosis <- factor(data_A$diagnosis, levels = c(0, 1), labels = c("No RA", "RA"))

# List of outcome variables
outcome_vars <- colnames(Y)

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
  df_plot <- data_A[, c("time", "diagnosis", outcome_var)]
  names(df_plot)[3] <- "outcome"  # Rename the outcome column for consistency
  
  # Fit smooth.spline for 'RA'
  df_ra <- df_plot[df_plot$diagnosis == "RA", ]
  spline_ra <- smooth.spline(df_ra$time, df_ra$outcome, df = 4)
  
  # Fit smooth.spline for 'no_RA'
  df_no_ra <- df_plot[df_plot$diagnosis == "No RA", ]
  spline_no_ra <- smooth.spline(df_no_ra$time, df_no_ra$outcome, df = 4)
  
  # Create a sequence of time points for smooth curves
  time_seq <- seq(min(data_A$time, na.rm = TRUE), max(data_A$time, na.rm = TRUE), length.out = 1000)
  
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
    labs(title = paste(outcome_var, "Serum Levels Relative to RA Diagnosis (Sample A)"),
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

ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
          plot_list[[4]], plot_list[[5]], plot_list[[6]], nrow = 3, ncol = 2, 
          legend = "right", align = "v", common.legend = TRUE)


### Sample B

source("~/Github/ra-biomarker/clean-data-B.R")

data_B$diagnosis <- factor(data_B$diagnosis, levels = c("Control", "Case"), labels = c("No RA", "RA"))
biomarkers <- c("aptivaccp3igg_≥5#00flu", "aptiva_acpafsiggvimentin2_≥5#00au",
                "aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
                "aptivaccp3iga_≥5#00flu", "aptiva_acpafsigavimentin2_≥5#00au",
                "aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")

biomarkers_labels <- c("anti-CCP3 (IgG)","anti-citVim2 (IgG)", "anti-citFib (IgG)","anti-citHis1 (IgG)",
                       "anti-CCP3 (IgA)","anti-citVim2 (IgA)","anti-citFib (IgA)","anti-citHis1 (IgA)")

outcome_colors <- brewer.pal(12, "Paired")

names(outcome_colors) <- biomarkers

# Define colors for 'RA' and 'no_RA'
group_colors <- c("RA" = "dodgerblue4", "No RA" = "grey")  # Blue for RA, gray for no_RA

# Transparency for scatter points
point_alpha <- 0.5

# Initialize a list to store individual plots
plot_list <- list()

# Loop through each outcome variable to create plots
for (i in seq_along(biomarkers)) {
  outcome_var <- biomarkers[i] 
  bio_label<-biomarkers_labels[i] # instead of full biomarker name use shortened version 
  color <- outcome_colors[outcome_var]
  
  # Prepare data for plotting
  df_plot <- data_B[, c("t_yrs", "diagnosis", outcome_var)]
  names(df_plot)[3] <- "outcome"  # Rename the outcome column for consistency
  #df_plot$outcome <- log(ifelse(df_plot[[outcome_var]] == 0, 0.00001, df_plot[[outcome_var]])) # log transform the outcome
  
  
  # Fit smooth.spline for 'RA'
  df_ra <- df_plot[df_plot$diagnosis == "RA", ]
  spline_ra <- smooth.spline(df_ra$t_yrs, df_ra$outcome, df = 4)
  
  # Fit smooth.spline for 'no_RA'
  df_no_ra <- df_plot[df_plot$diagnosis == "No RA", ]
  spline_no_ra <- smooth.spline(df_no_ra$t_yrs, df_no_ra$outcome, df = 4)
  
  # Create a sequence of time points for smooth curves
  time_seq <- seq(min(data_B$t_yrs, na.rm = TRUE), max(data_B$t_yrs, na.rm = TRUE), length.out = 1000)
  
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
               aes(x = t_yrs, y = outcome, color = diagnosis), 
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
    labs(
      x = "Time Prior to Diagnosis ",
      y = paste(bio_label)) +
    # Minimal theme for a clean look
    theme_minimal() +
    # Center the plot title and remove legend title
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # Add the plot to the list
  plot_list[[i]] <- p
  
}

ggarrange_args <- c(plot_list[1:8],  nrow = 2, ncol = 4, 
                    legend = "right", align = "v", common.legend = TRUE) 

combined_plot<-do.call(ggarrange, ggarrange_args)
final_plot <- annotate_figure(combined_plot, top = text_grob("Serum Levels Relative to RA Diagnosis (Sample B)", face = "bold", size = 14))
print(final_plot)
