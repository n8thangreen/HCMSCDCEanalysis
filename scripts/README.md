# R Analysis and Plotting Scripts

This directory contains the primary R scripts for data preparation, running the cost-effectiveness simulations, performing sensitivity analyses, and generating output figures.

## Files

### Data Preparation
* **`prep_study_data.R`**: Cleans the individual-level patient dataset, computes 5-year risk thresholds (>4% and >6%) using the Cox prediction model, and aggregates transition counts by strategy.
* **`raw-data-summary.R`**: Processes the raw dataset to output summaries of key events (follow-up duration, SCD, transplants, cardiovascular deaths).

### Cost-Effectiveness Simulations
* **`run_ce_analysis.R`**: **Unified runner script (Recommended)**. Allows running the complete cost-effectiveness analysis for either transvenous ICD or S-ICD, and under either a standard or tunnel-state model configuration, with simple parameter toggles at the top of the file. It runs simulations, executes BCEA, and saves CE planes and CEAC plots to the `images/` folder.
* **`ce_analysis_all_strat.R`**: Simulates CE outcomes across all decision rule stratifications over all deterministic scenarios from `scenario_input_values.csv` (used for generating Deterministic Sensitivity Analysis results).

### Sensitivity Analyses
* **`sensitivity_analysis_variance-based.R`**: Runs global variance-based sensitivity analysis (Sobol method) on the Markov model parameters.

### Plotting and Output Verification
* **`plot_pop_counts.R`**: Generates line plots showing state occupancy / population counts over time for each strategy.
* **`plot_posterior_hist.R`**: Generates histograms of the posterior transition probabilities from the JAGS output.
* **`plot_tornado_plots.R`**: Generates one-way Deterministic Sensitivity Analysis (DSA) tornado plots using Expected Incremental Benefit (EIB).
* **`check_ce_output_stats.R`**: Quick utility script to inspect summary statistics (means, standard deviations) of the cost and utility outputs.
