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

onedrive<- get_business_onedrive()
file_path <- "Attachments/KevinDat2.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
onedrive$download_file(
  src = file_path,
  dest = temp_file,
  overwrite = TRUE
)
dat_2raw<- read_xlsx(temp_file)
unlink(temp_file)

colnames(dat_2raw) <- tolower(colnames(dat_2raw))


all_biomarkers<-c("aptivaccp3iga_≥5#00flu","aptivaccp3igg_≥5#00flu","aptivapad1igg_≥5#00au","aptivapad2igg_≥5#00au",
              "aptivapad3igg_≥5#00au","aptivapad4igg_≥5#00au","aptivapad6igg_≥5#00au","aptiva_acpafsiggvimentin2_≥5#00au",
              "aptiva_acpafsigghistone2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au","aptiva_acpafsiggvimentin1_≥5#00au",
              "aptivapad1iga_≥5#00au","aptivapad2iga_≥5#00au","aptivapad3iga_≥5#00au","aptivapad4iga_≥5#00au",
              "aptivapad6iga_≥5#00au","aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigahistone2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au",
              "aptiva_acpafsigahistone1_≥5#00au","aptiva_acpafsigavimentin1_≥5#00au","quantalitecarpigg_≥20au",
              "quantalitecarpiga_≥20au","quantaflashrfiga_≥20#0cu_a","quantaflashrfigm_≥5#0iuml_a")

biomarkers<-c("aptivaccp3iga_≥5#00flu","aptivaccp3igg_≥5#00flu",
              "aptivapad1igg_≥5#00au","aptivapad4igg_≥5#00au",
              "aptiva_acpafsiggvimentin2_≥5#00au","aptiva_acpafsiggfibrinogen_≥5#00au","aptiva_acpafsigghistone1_≥5#00au",
              "aptivapad1iga_≥5#00au","aptivapad4iga_≥5#00au",
              "aptiva_acpafsigavimentin2_≥5#00au","aptiva_acpafsigafibrinogen_≥5#00au","aptiva_acpafsigahistone1_≥5#00au")
excu_biomarkers<-setdiff(all_biomarkers,biomarkers)
biomarkers<-excu_biomarkers # the excluded biomarkers are now the biomarkers 
dat_2 <- dat_2raw %>% 
  dplyr::rename(diagnosis=casecontrol, 
                study_id=masterstudyid,
                subj_id=masterstudyidnumeric,
                sampnum=sampordernum,
                age=ageserumsample,
                t_days=d_serum_ref) %>% 
  select(study_id, subj_id, sampnum,diagnosis, age,t_days,year_dref,
         gender,eversmoke,familyhxra,race_ethnic, all_of(biomarkers))

dat_2$subj_id <- ifelse(dat_2$diagnosis=="Control", paste0(dat_2$subj_id, "_1"), paste(dat_2$subj_id)) # control and case had the same subj ID
dat_2$subj_id <- as.character(as.integer(factor(dat_2$subj_id)))
dat_2$t_yrs <- dat_2$t_days/365 # changing from days to years

dat_2[,biomarkers] <- apply(dat_2[,biomarkers], 2, as.numeric)
dat_2 <- arrange(dat_2, as.numeric(subj_id), sampnum) %>% drop_na(all_of(biomarkers))
dat_2$diagnosis <- factor(dat_2$diagnosis, levels = c("Control", "Case"), labels = c("no_RA", "RA"))

outcome_colors <- c(brewer.pal(12, "Paired"),"#FF33FD","#33FFDE")
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
  color <- outcome_colors[outcome_var]
  
  # Prepare data for plotting
  df_plot <- dat_2[, c("t_yrs", "diagnosis", outcome_var)]
  names(df_plot)[3] <- "outcome"  # Rename the outcome column for consistency
  
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
    labs(title = "Serum Levels over Time with Smoothing Spline",
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


ggarrange_args <- c(plot_list[1:14],  nrow = 4, ncol = 4, 
                    legend = "right", align = "v", common.legend = TRUE) 
do.call(ggarrange, ggarrange_args)

