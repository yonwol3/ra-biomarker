library(corrplot)
library(tidyverse)

# Posterior correlation among the change-point (delta) parameters.

setwd("~/Documents/RA-Biomarker/")
model_tag <- "censcorr"

# ---- Sample A ----
load(sprintf("mcmc/mcmc_%s_A.RData", model_tag))
mcmc_A <- get(sprintf("mcmc_%s_A", model_tag))

biomarkers_A <- c("RF IgA", "RF IgM", "RF IgG",
                  "ACPA IgA", "ACPA IgM", "ACPA IgG")

mcmc_df_A <- as.data.frame(mcmc_A)
delta_cols_A <- paste0("delta[", 1:6, "]")
cor_delta_A  <- round(cor(mcmc_df_A[, delta_cols_A]), 2)
rownames(cor_delta_A) <- colnames(cor_delta_A) <- biomarkers_A

png(file = "figures/Sample A posterior correlation (delta).png",
    width = 10, height = 10, units = "in", res = 300)
corrplot(cor_delta_A, method = "color",
         title = "Cohort A: Posterior Correlation of Changepoints",
         mar = c(0, 0, 2, 0))
dev.off()

# ---- Sample B ----
load(sprintf("mcmc/mcmc_%s_B.RData", model_tag))
mcmc_B <- get(sprintf("mcmc_%s_B", model_tag))

biomarkers_B <- c("anti-CCP3 (IgG)", "anti-citVim2 (IgG)",
                  "anti-citFib (IgG)", "anti-citHis1 (IgG)",
                  "anti-CCP3 (IgA)", "anti-citVim2 (IgA)",
                  "anti-citFib (IgA)", "anti-citHis1 (IgA)")

mcmc_df_B <- as.data.frame(mcmc_B)
delta_cols_B <- paste0("delta[", 1:8, "]")
cor_delta_B  <- round(cor(mcmc_df_B[, delta_cols_B]), 2)
rownames(cor_delta_B) <- colnames(cor_delta_B) <- biomarkers_B

png(file = "figures/Sample B posterior correlation (delta).png",
    width = 10, height = 10, units = "in", res = 300)
corrplot(cor_delta_B, method = "color",
         title = "Cohort B: Posterior Correlation of Changepoints",
         mar = c(0, 0, 2, 0))
dev.off()
