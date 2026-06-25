# Bayesian Models (JAGS / BUGS)

This directory contains the WinBUGS/OpenBUGS/JAGS model specifications and R scripts used to estimate transition probabilities from discrete transition count data.

## Files

### Model Specifications (.txt)
* **`BUGS_code_model_all_strat.txt`**: Main BUGS model specification to estimate transition probabilities using multinomial likelihood across all risk stratifications.
* **`BUGS_code_model_all_strat_direct_dirichlet.txt`**: BUGS model utilizing a Dirichlet distribution directly as the prior for transition probabilities.
* **`BUGS_code_model_all_strat_rates.txt`**: BUGS model specified in terms of transition rates rather than probabilities directly.
* **`BUGS_code_model_with_no_icd.txt`**: BUGS model specification that includes transition paths for the 'no ICD' cohort.
* **`BUGS_code_all_strat_bivariate_logit.R`**: BUGS model for bivariate logit meta-analysis to model correlations between multiple parameters (e.g. shock and mortality).
* **`BUGS_code_best_model.txt`**: Combined best model specification using continuous Poisson transition rates with a bivariate normal prior to capture correlations between ICD and non-ICD death rates.

### R Runner Scripts (.R)
* **`BUGS_script_all_strat.R`**: Script to load aggregated data, run the main JAGS model via `R2jags`, check convergence, and save the posterior samples to `data/lambda.Rds`.
* **`BUGS_script_all_strat_direct_dirichlet.R`**: Script to run the direct Dirichlet conjugate model.
* **`BUGS_script_all_strat_bivariate_logit.R`**: Script to run the bivariate logit meta-analysis model.
* **`BUGS_script_with_no_icd.R`**: Script to run the model incorporating the 'no ICD' strategy and save cohort-specific probability arrays.
* **`BUGS_script_best_model.R`**: Runner script to run the combined best model using JAGS and save the posterior probability RDS files to `data/`.
