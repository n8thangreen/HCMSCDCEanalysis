# Cost-Effectiveness Modelling of S-ICD vs. TV-ICD in Hypertrophic Cardiomyopathy

## 1. Health Benefit Dynamics in Unselected HCM Populations

In general Hypertrophic Cardiomyopathy (HCM) populations, the absolute Quality-Adjusted Life Year (QALY) benefit of implanting an Implantable Cardioverter-Defibrillator (ICD) is mathematically small due to three primary drivers:

* **Low Baseline Risk:** The annual incidence of Sudden Cardiac Death (SCD) in unselected cohorts is historically around 1%. Without filtering for high-risk patients, the number needed to treat to save one life is exceptionally high.
* **Unresolved Morbidity:** ICDs prevent fatal ventricular arrhythmias but do not treat the underlying cardiomyopathy, LVOT obstruction, or heart failure symptoms that primarily degrade quality of life.
* **Complication Disutilities:** HCM patients are often diagnosed young, exposing them to a lifetime risk of device complications (lead fractures, inappropriate shocks, systemic infections). The QALY penalties applied to these adverse events heavily offset the survival gains.

Cost-effectiveness relies heavily on risk stratification. Subcutaneous ICDs (S-ICD) and Transvenous ICDs (TV-ICD) also differ structurally; S-ICDs avoid vascular lead complications but lack Anti-Tachycardia Pacing (ATP), altering the distribution of complication-related disutilities.

## 2. Device Longevity and Replacement Distributions

Modelling device replacements at fixed 8- or 10-year timepoints forces an artificial step-function into the Markov trace, causing unnatural, simultaneous spikes in costs and peri-operative mortality. Best practice handles this via:

* **Structural Time-Dependency:** Parametric survival curves (typically Weibull) are used to capture the monotonically increasing hazard of battery depletion, smoothing replacement events over time.
* **Parameter Uncertainty:** During Probabilistic Sensitivity Analysis (PSA), uncertainty is propagated by sampling the survival parameters directly.

## 3. Implementing Markov Tunnel States

Tunnel states bypass the standard Markov memoryless assumption by indexing transition probabilities to the time since the device was implanted.

1. **Tunnel Sequence:** The primary "ICD Active" state is expanded into a series of temporary tunnel states equal to the maximum possible lifespan of the device (`ICD_cycle_1`, `ICD_cycle_2` ... `ICD_cycle_max`).
2. **Cycle-Specific Probabilities:** For each tunnel state $t$, the probability of replacement is derived from the discrete-time hazard using the survival function $S(t)$:

$$P(\text{replace}_{t}) = 1 - \frac{S(t)}{S(t-1)}$$


3. **Transition Routing:** Patients transition to the `Replacement` state based on $P(\text{replace}_{t})$. Surviving patients who receive a replacement transition directly back to `ICD_cycle_1` in the subsequent cycle, resetting the time-to-event clock.

## 4. R Implementation: Base Replacement Probabilities

The following code calculates the cycle-dependent replacement probability vectors for tunnel states using a Weibull distribution parametrised from median device lifespans.

```r
# Calculate Weibull scale parameter from median
calc_weibull_scale <- function(median_time, shape) {
  median_time / (log(2)^(1/shape))
}

shape_param <- 5
tv_scale <- calc_weibull_scale(median_time = 10.0, shape = shape_param)
sicd_scale <- calc_weibull_scale(median_time = 8.7, shape = shape_param)

# Calculate conditional transition probabilities for tunnel states
calc_replacement_probs <- function(max_cycles, shape, scale) {
  cycles <- 1:max_cycles
  surv_t <- pweibull(cycles, shape = shape, scale = scale, lower.tail = FALSE)
  surv_t_minus_1 <- pweibull(cycles - 1, shape = shape, scale = scale, lower.tail = FALSE)
  
  # P(event in cycle t | survived to cycle t-1)
  ifelse(surv_t_minus_1 == 0, 1, 1 - (surv_t / surv_t_minus_1))
}

max_years <- 15
p_tv_replace <- calc_replacement_probs(max_years, shape_param, tv_scale)
p_sicd_replace <- calc_replacement_probs(max_years, shape_param, sicd_scale)

```

## 5. Dynamic Parameter Sampling for PSA

During PSA, Weibull parameters must be sampled dynamically per iteration. The optimal method depends on data availability.

### Joint Sampling from Fitted Data (Gold Standard)

If patient-level data is available, parameters are drawn from a multivariate normal distribution on the log scale to preserve covariance.

```r
library(MASS)
library(flexsurv)

# Fit model and extract variance-covariance matrix
fit_sicd <- flexsurvreg(Surv(time, status) ~ 1, data = device_data, dist = "weibull")
log_estimates <- fit_sicd$res[, "est"]
vcov_matrix <- fit_sicd$vcov

# Draw joint samples for PSA
n_sims <- 1000
psa_log_samples <- mvrnorm(n = n_sims, mu = log_estimates, Sigma = vcov_matrix)
psa_shapes <- exp(psa_log_samples[, "shape"])
psa_scales <- exp(psa_log_samples[, "scale"])

```

### Sampling from Published Medians (Pragmatic Approach)

If relying on published summary statistics, a distribution is placed over the median itself, and the scale parameter is dynamically recalculated per iteration alongside a sampled shape parameter.

```r
n_sims <- 1000

# Sample the median (e.g., S-ICD median 8.7 years, 95% CI: 7.9 to 9.5)
median_mean <- 8.7
median_se <- (9.5 - 7.9) / (2 * 1.96) 
psa_medians <- rnorm(n_sims, mean = median_mean, sd = median_se)

# Sample the shape parameter using a Gamma distribution
shape_mean <- 5
shape_se <- 0.5
psa_shapes <- rgamma(n_sims, 
                     shape = (shape_mean / shape_se)^2, 
                     rate = shape_mean / (shape_se^2))

# Calculate tunnel probabilities per iteration
max_years <- 15
psa_prob_matrix <- matrix(NA, nrow = n_sims, ncol = max_years)

for (i in 1:n_sims) {
  scale_i <- calc_weibull_scale(median_time = psa_medians[i], shape = psa_shapes[i])
  psa_prob_matrix[i, ] <- calc_replacement_probs(max_years, psa_shapes[i], scale_i)
}

```