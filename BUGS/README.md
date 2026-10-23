# Bayesian Models (JAGS / BUGS)

This directory contains the JAGS/BUGS model specifications and R runner scripts used to estimate transition probabilities from cohort count data.

## Files

### Model Specifications (.txt)
* **`model_multinomial_dirichlet.txt`**: Uses a **Multinomial** likelihood (`dmulti`) with a **Dirichlet** prior (`ddirch`) to estimate independent transition probabilities for each cohort (strategy).
* **`model_multinomial_dirichlet_noicd.txt`**: Same as above, but explicitly models the additional "no ICD" strategy/cohort.
* **`model_direct_dirichlet_conjugate.txt`**: Samples directly from the conjugate **Dirichlet posterior** distribution based on prior and transition counts, bypassing the need for an explicit JAGS likelihood statement (`dmulti`).
* **`model_poisson_rates.txt`**: Uses a **Poisson** likelihood on transition counts to estimate continuous transition rates, which are then converted to discrete cycle probabilities using a competing-risks formula.
* **`model_binomial_bivariate_logit_death.txt`**: The **best-performing model**. It splits transitions into independent **Binomial** likelihoods (`dbin`) to allow a **bivariate normal prior** on logit-probabilities of death (`p_icd_death` and `p_noicd_death`), capturing their correlation across risk strata.

### R Runner Scripts (.R)
* **`run_multinomial_dirichlet.R`**: Runs `model_multinomial_dirichlet.txt` and saves output to `data/lambda.Rds`.
* **`run_multinomial_dirichlet_noicd.R`**: Runs `model_multinomial_dirichlet_noicd.txt` and saves strategy-specific outputs to `data/`.
* **`run_direct_dirichlet_conjugate.R`**: Runs `model_direct_dirichlet_conjugate.txt` (the direct conjugate updater).
* **`run_binomial_bivariate_logit_death.R`**: Runs the best model `model_binomial_bivariate_logit_death.txt` and saves outputs to `data/`.
