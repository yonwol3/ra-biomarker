library(corrplot)
library(tidyverse)

# ---- Sample A: Posterior correlation from MCMC draws ----

setwd("~/Documents/RA-Biomarker/")
load("mcmc/mcmc_trunc_A.RData")

# Biomarker labels corresponding to gamma[1:6] / delta[1:6]
biomarkers_A <- c("RF IgA", "RF IgM", "RF IgG",
                  "ACPA IgA", "ACPA IgM", "ACPA IgG")

# Convert MCMC matrix to data frame
mcmc_df_A <- as.data.frame(mcmc_trunc_A)

# Posterior correlation among delta parameters
delta_cols_A <- paste0("delta[", 1:6, "]")
cor_delta_A  <- round(cor(mcmc_df_A[, delta_cols_A]), 2)
rownames(cor_delta_A) <- colnames(cor_delta_A) <- biomarkers_A

png(file = "figures/Sample A posterior correlation (delta).png",
    width = 10, height = 10, units = "in", res = 300)
corrplot(cor_delta_A, method = "color",
         title = "Sample A: Posterior Correlation of Changepoints",
         mar = c(0, 0, 2, 0))
dev.off()


# ---- Sample B: Posterior correlation from MCMC draws ----

setwd("~/Documents/RA-Biomarker/")
load("mcmc/mcmc_trunc_B.RData")

biomarkers_B <- c("anti-CCP3 (IgG)", "anti-citVim2 (IgG)",
                  "anti-citFib (IgG)", "anti-citHis1 (IgG)",
                  "anti-CCP3 (IgA)", "anti-citVim2 (IgA)",
                  "anti-citFib (IgA)", "anti-citHis1 (IgA)")

mcmc_df_B <- as.data.frame(mcmc_trunc_B)

# Posterior correlation among delta parameters
delta_cols_B <- paste0("delta[", 1:8, "]")
cor_delta_B  <- round(cor(mcmc_df_B[, delta_cols_B]), 2)
rownames(cor_delta_B) <- colnames(cor_delta_B) <- biomarkers_B

png(file = "figures/Sample B posterior correlation (delta).png",
    width = 10, height = 10, units = "in", res = 300)
corrplot(cor_delta_B, method = "color",
         title = "Sample B: Posterior Correlation of Changepoints",
         mar = c(0, 0, 2, 0))
dev.off()
