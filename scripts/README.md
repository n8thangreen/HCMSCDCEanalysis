# R Analysis and Plotting Scripts

This directory contains the primary R scripts for data preparation, running the cost-effectiveness simulations, performing sensitivity analyses, and generating output figures.

## Files

### Data Preparation
* **`prep_study_data.R`**: Cleans the individual-level patient dataset, computes 5-year risk thresholds (>4% and >6%) using the Cox prediction model, and aggregates transition counts by strategy.
* **`raw-data-summary.R`**: Processes the raw dataset to output summaries of key events (follow-up duration, SCD, transplants, cardiovascular deaths).

### Cost-Effectiveness Simulations
* **`run_ce_analysis.R`**: **Unified runner script (Recommended)**. Allows running the complete cost-effectiveness analysis for either transvenous ICD or S-ICD, and under either a standard or tunnel-state model configuration, with simple parameter toggles at the top of the file. It runs simulations, executes BCEA, and saves CE planes and CEAC plots to the `images/` folder.
* **`run_combined_device_psa.R`**: Runs a combined TV-ICD and S-ICD probabilistic sensitivity analysis (PSA) comparing both devices across their respective tunnel-state models (under baseline parameters) on the same CE Plane and CEAC plots.
* **`run_baseline_model_comparison.R`**: Runs the 4 main combinations of device types and model configurations (TV-ICD standard/tunnel, S-ICD standard/tunnel) using baseline point estimate parameters and outputs a summary table to `data/all_ce_model_results.csv`.
* **`ce_analysis_all_strat_s-icd.R`**: Simulates CE outcomes across all decision rule stratifications over all deterministic scenarios from `scenario_input_values.csv` (used for generating Deterministic Sensitivity Analysis results).

### Paper Replication
* **`run_original_icd_analysis.R`**: Replicates the deterministic base-case cost-effectiveness results for TV-ICD (horizon = 12 years) as reported in Table 2 of the journal paper.
* **`run_original_icd_analysis_psa.R`**: Replicates the probabilistic sensitivity analysis (PSA) results for TV-ICD (horizon = 12 years) as reported in Table 2 and Figure 2 of the paper.

### Sensitivity Analyses
* **`sensitivity_analysis_variance-based.R`**: Runs global variance-based sensitivity analysis (Sobol method) on the Markov model parameters.

### Plotting and Output Verification
* **`plot_pop_counts.R`**: Generates line plots showing state occupancy / population counts over time for each strategy.
* **`plot_posterior_hist.R`**: Generates histograms of the posterior transition probabilities from the JAGS output.
* **`plot_tornado_plots.R`**: Generates one-way Deterministic Sensitivity Analysis (DSA) tornado plots using Expected Incremental Benefit (EIB).
