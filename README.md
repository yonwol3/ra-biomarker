# A Bayesian Multivariate Segmented Regression Analysis for Biomarker Detection in Rheumatoid Arthritis

The repository contains code, figures, and tables applied to 'Detecting Change-Points in Preclinical Rheumatoid Arthritis Biomarkers using Bayesian Multivariate Segmented Regression.'

## Contents

### [`figures`](https://github.com/yonwol3/ra-biomarker/tree/main/figures) Folder
- Contains all the figures generated for our analysis

### [`stan`](https://github.com/yonwol3/ra-biomarker/tree/main/stan) Folder
- Contains STAN Models (for both truncation, censored, and binarized models)

### [`simulation`](https://github.com/yonwol3/ra-biomarker/tree/main/simulation) Folder
- Contains code for simulating exposure-outcome relationships. Purely pedagological, but useful for testing methods.

### R Scripts 
- [`clean-data-A.R`](https://github.com/yonwol3/ra-biomarker/blob/main/clean-data-A.R): data cleaning for sample A.
- [`clean-data-B.R`](https://github.com/yonwol3/ra-biomarker/blob/main/clean-data-B.R): data cleaning for sample B.
- [`mcmc-stan.R`](https://github.com/yonwol3/ra-biomarker/blob/main/mcmc-stan.R): Interface with STAN to draw MCMC samples for all models.
- [`hpd.R`](https://github.com/yonwol3/ra-biomarker/blob/main/hpd.R): Function for constructing highest posterior density credible intervals.
- [`loess-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/loess-plots.R): Code to generate the loess plots.
- [`tables.R`](https://github.com/yonwol3/ra-biomarker/blob/main/tables.R): includes  Code used to generate Table 1.
- [`truncated-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/truncated-plots.R): code to generate the posterior density plots and summaries for the truncated model (both sample A and B biomarkers) 
- [`censored-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/censored-plots.R): code to generate the posterior density plots and summaries for the right-censored model (both sample A and B biomarkers) 
- [`binary-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/binary-plots.R): code to generate the posterior density plots and summaries after binarizing outcome model (both sample A and B biomarkers) 

## References