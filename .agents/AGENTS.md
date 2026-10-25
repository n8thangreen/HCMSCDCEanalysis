# Project-Scoped Rules & Context: HCMSCDCEanalysis

This file provides critical context, model structures, and implementation rules for the **HCMSCDCEanalysis** R package. Future agents working in this workspace should adhere to these instructions.

---

## 1. Project Domain & Context
This package simulates and analyzes the cost-effectiveness of **Transvenous ICD (TV-ICD / VTICD)** vs. **Subcutaneous ICD (S-ICD / SICD)** strategies in hypertrophic cardiomyopathy (HCM) patients using Markov models.

### Models:
1. **Standard 6-State Model**: 
   * Compares cohorts across 6 health states.
   * **Rule**: Requires `inp$t_repl_max` to be strictly `1` to avoid dimension mismatch errors between population matrices and state cost/utility vectors.
2. **Tunnel State Model (Device Replacement Cycles)**:
   * Expands the model to track time-since-implant (tunneling).
   * **State space size**: `2 * t_repl_max + 4` (where `t_repl_max` is the number of tunnel states representing battery lifespan years).

---

## 2. Core Package Rules & Best Practices

### A. MCMC Output & Lazy-loading Workaround
* **Problem**: Both `data/jags_output.RData` and `data/jags_output_all_strat_direct.RData` define a JAGS result object named `mm1`. Since R packages compile datasets into a single lazy-loaded database, one `mm1` object silently overwrites the other.
* **Rule**: **NEVER** use `data("jags_output")` in workspace scripts. Instead, load the file directly using:
  ```R
  load(here::here("data/jags_output.RData"))
  ```

### B. Function Signatures & `safe_call` Compatibility
* **Background**: The simulation engine in `ce_sim` is designed to be backwards-compatible with older scripts/tests that call cost or utility functions with fewer arguments.
* **Rule**: Always use the `safe_call` wrapper in `ce_sim.R` to invoke cost and utility sub-functions (`rc_init`, `rc_unit`, `re_init`, `re_unit`), as it automatically maps arguments according to the function's signature.

### C. JAGS/BUGS Parameter Initialization
* **Problem**: JAGS will fail with a `"Dimension mismatch in values supplied for rate"` if you attempt to initialize deterministic/unused entries in the `rate` array.
* **Rule**: When running the best model JAGS MCMC script (`BUGS/run_binomial_bivariate_logit_death.R`), define initial values for `rate` using a dimension of `c(4, 4, 6)` and fill only the active stochastic indices (`[i, 1, 2]` and `[i, 4, 5]`), leaving others as `NA`.

### D. Parameter csv Files (`inst/extdata`)
* The csv files `scenario_input_values.csv`, `input_distns.csv`, and `input_values.csv` contain a `label` or `Scenario_Label` column with values `tvicd` (old ICD) and `sicd` (new ICD) to hold both model configurations concurrently.
* **Rule**: Preserve these columns and use them for filtering or indexing scenarios to ensure compatibility with both device types.
