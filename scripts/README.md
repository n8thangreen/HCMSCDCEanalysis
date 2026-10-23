# R Analysis and Plotting Scripts (30-Year Tunnel-State S-ICD vs. TV-ICD Analysis)

This directory contains the active R scripts for data preparation, running the 30-year tunnel-state Markov cost-effectiveness simulations (comparing **TV-ICD** with a 10-year replacement cycle vs. **S-ICD** with an 8-year replacement cycle), performing sensitivity analyses, and generating manuscript figures.

Legacy replication scripts from the 2024 *EHJ-QCCO* 12-year non-tunnel analysis have been moved to [`scripts/archive/`](archive/).

---

## Active Pipeline

### 1. Data Preparation
* **`prep_study_data.R`**: Cleans the individual-level patient dataset (`data_set1`), computes 5-year SCD risk thresholds (`>4%` and `>6%`) and risk factor counts, and saves `data/ipd_risk.RData` and `data/transdat.RData`.

### 2. Primary Cost-Effectiveness Analyses (TV-ICD vs. S-ICD)
* **`run_baseline_model_comparison.R`**: Runs the deterministic base-case 30-year tunnel-state model for both **TV-ICD** (10-year replacement cycle) and **S-ICD** (8-year replacement cycle) across the 4 core strategies (`Risk > 4%`, `Risk > 6%`, `Observed`, `No ICD`) and saves `data/all_ce_model_results.csv` (Table 2 in the manuscript).
* **`run_combined_device_psa.R`**: Runs the joint 6-strategy Probabilistic Sensitivity Analysis (PSA) comparing TV-ICD and S-ICD across `Cox > 4%`, `Cox > 6%`, and `Observed` against `TV-ICD: Observed` as reference, generating `images/ce_plane_combined.png` and `images/ceac_combined.png`.
* **`run_all_strat_comparison.R`** & **`plot_ce_all_strat.R`**: Runs the head-to-head S-ICD vs. TV-ICD comparison across all 8 risk stratification rules (`data/ce_analysis_all_strat_s-icd_results.RData`, Table 3 in the manuscript) and generates the INMB forest plot (`images/inmb_strategies_forest.png`) and strategy CE plane (`images/ce_plane_strategies.png`).

### 3. Sensitivity Analyses
* **`run_ce_analysis_dsa.R`** & **`plot_tornado_plots.R`**: Runs the one-way Deterministic Sensitivity Analysis (DSA) across all 36 parameter scenarios in `inst/extdata/scenario_input_values.csv` (`data/res_dsa.RData`) and generates Expected Incremental Benefit (EIB) tornado plots for both TV-ICD (`images/tornado_plots_EIB.png`) and S-ICD (`images/tornado_plots_EIB_sicd.png`).
* **`sensitivity_analysis_variance-based.R`**: Runs global variance-based sensitivity analysis (Sobol first-order indices) and generates `images/var_explain_barplot.png`.

### 4. State Occupancy & Bayesian Model Diagnostics
* **`plot_pop_counts.R`**: Generates 30-year cohort state occupancy plots for TV-ICD (`images/state_pop_over_time_tvicd_tunnel.png`) and S-ICD (`images/state_pop_over_time_sicd_tunnel.png`).
* **`plot_posterior_hist.R`** & **`plot_posterior_checking.R`**: Generates posterior transition probability density plots (`images/post_hist.png`) and posterior predictive check plots (`images/posterior_predictive_checks.png`).
