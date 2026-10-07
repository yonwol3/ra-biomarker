# A Bayesian Multivariate Segmented Regression Analysis for Biomarker Detection in Rheumatoid Arthritis

The repository contains code, figures, and tables applied to 'Detecting Change-Points in Preclinical Rheumatoid Arthritis Biomarkers using Bayesian Multivariate Segmented Regression.'

## Models

Both models are fit on the log scale from the same data list (`standata_cc_*`) with the same priors, so the primary and the sensitivity differ **only** in the residual correlation structure.

- **Primary:** `stan/change-point-cens-corr.stan` — values above the assay upper limit of detection are treated as right-censored via data augmentation (latent exceedances bounded below by the detection limit), with an unstructured residual correlation (LKJ prior on the Cholesky factor). This keeps censoring and inter-biomarker correlation in the model simultaneously.
- **Sensitivity:** `stan/change-point-cens.stan` — same censored likelihood, evaluated with the `1 - Phi(z)` survival term and a `log(fmax(pi, 1e-10))` guard, with a **diagonal** residual covariance.
- **Retired:** `stan/change-point-trunc.stan` — the truncated-likelihood model. The file is kept for provenance; it is commented out in `mcmc-stan.R` and is no longer fit or reported.

Sample B values reported as exactly 0 are set to half the minimum positive value for that biomarker in `clean-data-B.R` before the log transform.

## Contents

### [`figures`](https://github.com/yonwol3/ra-biomarker/tree/main/figures) Folder
- Contains all the figures generated for our analysis

### [`stan`](https://github.com/yonwol3/ra-biomarker/tree/main/stan) Folder
- Contains the STAN models (censored + correlated primary, diagonal censored sensitivity, retired truncated model, and the simulation models)

### [`simulation`](https://github.com/yonwol3/ra-biomarker/tree/main/simulation) Folder
- Contains code for simulating exposure-outcome relationships. Purely pedagological, but useful for testing methods.

### R Scripts
- [`clean-data-A.R`](https://github.com/yonwol3/ra-biomarker/blob/main/clean-data-A.R): data cleaning for sample A.
- [`clean-data-B.R`](https://github.com/yonwol3/ra-biomarker/blob/main/clean-data-B.R): data cleaning for sample B.
- [`mcmc-stan.R`](https://github.com/yonwol3/ra-biomarker/blob/main/mcmc-stan.R): Interface with STAN to draw MCMC samples for the primary (censored + correlated) and sensitivity (diagonal censored) models. Writes `mcmc/mcmc_censcorr_{A,B}.RData` and `mcmc/mcmc_cens_{A,B}.RData`.
- [`hpd.R`](https://github.com/yonwol3/ra-biomarker/blob/main/hpd.R): Function for constructing highest posterior density credible intervals.
- [`loess-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/loess-plots.R): Code to generate the smoothing-spline plots of log serum levels against time to diagnosis (Figures 1 and 2).
- [`tables.R`](https://github.com/yonwol3/ra-biomarker/blob/main/tables.R): includes  Code used to generate Table 1.
- [`censored-corr-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/censored-corr-plots.R): **primary** results — posterior change-point densities (`figures/censcorr_change-point-dens_{A,B}.png`), composite-time summary tables (`tables/censcorr_summary_{A,B}.csv`), and the combined main-text Figure 3 panel.
- [`censored-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/censored-plots.R): **sensitivity** results — the same densities and summaries for the diagonal censored model (`figures/cens_change-point-dens_{A,B}.png`, `tables/censored_summary_{A,B}.csv`). Does not touch Figure 3.
- [`correlation-plots.R`](https://github.com/yonwol3/ra-biomarker/blob/main/correlation-plots.R): posterior change-point (delta) correlation matrices (Supplement), parameterised by `model_tag` (default "censcorr").

Each script re-sources the cleaning scripts and leaves objects named after base functions (`mean`, `t`, `gamma`) in the global environment. **Run each script in a fresh R session**, in the order: `mcmc-stan.R` → plotting/table scripts.

## References

https://www.medrxiv.org/content/10.64898/2026.05.22.26353892v2
