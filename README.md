A Bayesian Multivariate Segmented Regression Analysis for Biomarker Detection in Rheumatoid Arthritis
=====================================================================================================

The repository contains code, figures, and tables applied to 'Detecting Change-Points in Preclinical Rheumatoid Arthritis Biomarkers using Bayesian Multivariate Segmented Regression.'

## Contents

### [`figures`]() Folder
- contains all the figures generated for our analysis

### [`stan`]() Folder
- contains stan code (for both truncation, censored, and binarized models)

### R Scripts 
- [`clean-data-A.R`](): data cleaning for sample A.
- [`clean-data-B.R`](): data cleaning for sample B.
- [`mcmc-stan.R`](): Interface with STAN to draw MCMC samples for all models.
- [`hpd.R`](): Function for constructing highest posterior density credible intervals.
- [`loess-plots.R`](): Code to generate the loess plots.
- [`tables.R`](): includes  Code used to generate Table 1.
- [`truncated-plots.R`](): code to generate the posterior density plots and summaries for the truncated model (both sample A and B biomarkers) 
- [`censored-plots.R`](): code to generate the posterior density plots and summaries for the right-censored model (both sample A and B biomarkers) 
- [`binary-plots.R`](): code to generate the posterior density plots and summaries after binarizing outcome model (both sample A and B biomarkers) 

## References
