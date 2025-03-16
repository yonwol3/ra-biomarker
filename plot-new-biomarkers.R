library(tidyverse)
library(plyr)
library(splines)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

#--------------------------------#
# Loess plot of new 12 biomarkers
#--------------------------------#

source("clean-data-new.R")
dat_2$diagnosis <- factor(dat_2$diagnosis, levels = c("Control", "Case"), labels = c("no_RA", "RA"))

biomarkers<-c("aptivaccp3iga_≥5#00flu","aptivaccp3igg_≥5#00flu",
                          "aptivapad1igg_≥5#00au","aptivapad4igg_≥5#00au",
                          "aptiva_acpafsiggvimentin2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
                          "aptivapad1iga_≥5#00au","aptivapad4iga_≥5#00au",
                          "aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
biomarkers_labels<- sub("^aptiva", "", biomarkers)
biomarkers_labels<- sub("_≥5#00au", "", biomarkers_labels)
biomarkers_labels<-sub("_≥5#00flu","",biomarkers_labels)
biomarkers_labels<-sub("_","",biomarkers_labels)

outcome_colors <- brewer.pal(12, "Paired")
names(outcome_colors) <- biomarkers

# Define colors for 'RA' and 'no_RA'
group_colors <- c("RA" = "dodgerblue4", "no_RA" = "grey")  # Blue for RA, gray for no_RA

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
  df_plot <- dat_2[, c("t_yrs", "diagnosis", outcome_var)]
  names(df_plot)[3] <- "outcome"  # Rename the outcome column for consistency
  #df_plot$outcome <- log(ifelse(df_plot[[outcome_var]] == 0, 0.00001, df_plot[[outcome_var]])) # log transform the outcome

  
  # Fit smooth.spline for 'RA'
  df_ra <- df_plot[df_plot$diagnosis == "RA", ]
  spline_ra <- smooth.spline(df_ra$t_yrs, df_ra$outcome, df = 4)
  
  # Fit smooth.spline for 'no_RA'
  df_no_ra <- df_plot[df_plot$diagnosis == "no_RA", ]
  spline_no_ra <- smooth.spline(df_no_ra$t_yrs, df_no_ra$outcome, df = 4)
  
  # Create a sequence of time points for smooth curves
  time_seq <- seq(min(dat_2$t_yrs, na.rm = TRUE), max(dat_2$t_yrs, na.rm = TRUE), length.out = 1000)
  
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


ggarrange_args <- c(plot_list[1:12],  nrow = 3, ncol = 4, 
                    legend = "right", align = "v", common.legend = TRUE
                    ) 
combined_plot<-do.call(ggarrange, ggarrange_args)
final_plot <- annotate_figure(combined_plot, 
                              top = text_grob("Serum Levels over Time with Smoothing Spline (Sample 2)", 
                                              face = "bold", size = 14))
print(final_plot)

#-------------------------------------------#  
# Density plots from posterior distribution of
# model parameters (for 12 biomarkers)
#--------------------------------------------#

onedrive<- get_business_onedrive()
file_path <- "Attachments/mcmc_new_trunc.RData"
temp_file <- tempfile(fileext = ".RData")
onedrive$download_file(
  src = file_path,
  dest = temp_file,
  overwrite = TRUE
)
mcmc_new_trunc<- get(load(temp_file)[1])
mcmc<-mcmc_new_trunc[[1]]
kappa<-mcmc[,13:24] # change points only (no gamma (rate of change))

png("figures/change-point-dens-new-biomarkers.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(kappa[,1]), lwd = 2,
     col = outcome_colors[1], ylab = "Posterior Density", xlab = "Years Prior to Diagnosis",
     ylim = c(0, 0.8),
     xlim = c(-20, 5),
     main = "Change Point Densities (Sample 2)")

for (i in 2:12) {
  lines(density(kappa[,i]), lwd = 2, col = outcome_colors[i])
}

abline(v = 0, lty = 2, col = "blue")
abline(h = 0, lty = 1, col = "black")
grid()

legend("topleft", 
       legend = biomarkers_labels,
       col = outcome_colors, 
       lwd = rep(2, 12),
       cex = 1)

dev.off()


#-------------------------------------------#
# calculate mean and 95% Credible Interval 
#-------------------------------------------#
credible_res <- data.frame(biomarker=character(0),mean = numeric(0), 
                           credible_lower = numeric(0),
                           credible_upper = numeric(0))

for (i in 1:ncol(kappa)) {
  credible_res[i, "mean"] <- round(mean(kappa[, i]), 2)         # mean 
  credible_res[i, "credible_lower"] <- round(quantile(kappa[, i], 0.025), 2)  # 2.5th quantile
  credible_res[i, "credible_upper"] <- round(quantile(kappa[, i], 0.975), 2)  # 97.5th quantile
  credible_res[i, "biomarker"]<- biomarkers_labels[i]
}
credible_res<-arrange(credible_res,mean)
credible_res<-credible_res %>% 
  mutate(`95% Credible Interval`=paste("[", credible_lower, ", ", credible_upper, "]"))%>%
  select(biomarker, mean, `95% Credible Interval`)

write.csv(credible_res,"../../new_cp_ci.csv")

