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

# Reading data cleaning R file
source("~/Github/ra-biomarker/clean-dat.R")

# Set working directory
setwd("~/Dropbox/Projects/RA-Biomarker/")

# Naming serum antibodies correctly
colnames(Y) <- c("RF IgA", "RF IgM", "RF IgG", "ACPA IgA", "ACPA IgM", "ACPA IgG")

## Descriptive Spline Plot

png(filename = "images/spline-plots-horiz.png",
    width = 12,
    height = 8,
    units = "in",
    res = 300)

par(mfrow = c(2,3))

# smoothing splines - igarfconc_
sm_spline1 <- smooth.spline(raDat$t_yrs, logY[,1], df = 4)
plot(raDat$t_yrs, logY[,1], col="grey",
     xlab = "Time before Diagnosis", ylab = "RF-IgA",
     main = "RF-IgA Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline1, col = "darkolivegreen4", lwd = 2)

#smoothing splines - igmrfconc_
sm_spline2 <- smooth.spline(raDat$t_yrs, logY[,2], df = 4)
plot(raDat$t_yrs, logY[,2], col = "grey",
     xlab = "Time before Diagnosis", ylab = "RF-IgM", 
     main = "RF-IgM Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline2, col = "darkorange4", lwd = 2)

#smoothing splines - iggrfconc_
sm_spline3 <- smooth.spline(raDat$t_yrs, logY[,3], df = 4)
plot(raDat$t_yrs, logY[,3], col = "grey",
     xlab = "Time before Diagnosis", ylab="RF-IgG", 
     main = "RF-IgG Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline3, col= "firebrick4", lwd = 2)

#smoothing splines - igaccpavgconc
sm_spline4 <- smooth.spline(raDat$t_yrs, logY[,4], df = 4)
plot(raDat$t_yrs, logY[,4],col = "grey",
     xlab = "Time before Diagnosis", ylab="ACPA-IgA",
     main="ACPA-IgA Serum Levels over Time with \nSmoothing Spline")
abline(v=0, lty=3, col = "black")
lines(sm_spline4, col = "dodgerblue4", lwd = 2)

#smoothing splines - igmccpavgconc
sm_spline5 <- smooth.spline(raDat$t_yrs, logY[,5], df = 4)
plot(raDat$t_yrs, logY[,5], col="grey",
     xlab = "Time before Diagnosis", ylab = "ACPA-IgM",
     main = "ACPA-IgM Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline5, col = "mediumorchid4", lwd = 2)

sm_spline6 <- smooth.spline(raDat$t_yrs, logY[,6], df = 4)
plot(raDat$t_yrs, logY[,6], col="grey",
     xlab = "Time before Diagnosis", ylab = "ACPA-IgG",
     main = "ACPA-IgG Serum Levels over Time with \nSmoothing Spline")
abline(v = 0, lty = 3, col = "black")
lines(sm_spline6, col= "khaki4", lwd = 2)
dev.off()

## Changepoint Density Plot

load("mcmc/mcmc_b.RData")

kappa.names <- c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]", "kappa[5]", "kappa[6]")

kappa <- mcmc_b[[1]][,kappa.names]
colnames(kappa) <- colnames(Y)

png("images/change-point-dist.png", 
    width = 1000, 
    height = 1000,
    res = 100, 
    units = "px")

plot(density(kappa[,1]), lwd = 2,
     col = "darkolivegreen4", ylab = "Posterior Density", xlab = "kappa",
     ylim = c(0, 1.2),
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
abline(v = 0, lty = 3, col = "black")
legend("topleft", 
       legend = colnames(Y),
       col = c("darkolivegreen4", "darkorange4", "firebrick4", 
               "dodgerblue4", "mediumorchid4", "khaki4"), 
       lwd = c(2, 2, 2, 2, 2, 2),
       cex = 1)

dev.off()

## Diagnostics of model (a)

pdf("images/trace_a.pdf")
plot(mcmc_b)
dev.off()

mcmc <- as.matrix(mcmc[[1]])

png("images/acf_a.png", 
    width = 7000, 
    height = 6000,
    res = 100, 
    units = "px")
par(mfrow = c(4, 6), mar = c(1,1,1,1))
for (i in 1:ncol(mcmc))
  acf(mcmc[,i], ylab = colnames(mcmc)[i])
title("ACF Plots", outer = TRUE)
dev.off()

## Diagnostics of model (b)

pdf("images/trace_b.pdf")
plot(mcmc)
dev.off()

mcmc <- as.matrix(mcmc[[1]])

png("images/acf_b.png", 
    width = 7000, 
    height = 6000,
    res = 100, 
    units = "px")
par(mfrow = c(4, 6), mar = c(1,1,1,1))
for (i in 1:ncol(mcmc))
  acf(mcmc[,i], ylab = colnames(mcmc)[i])
title("ACF Plots", outer = TRUE)
dev.off()
