library(R2jags)
library(here)

# 1. Load Data
load("data/transdat.RData")

n_intervs <- length(transdat)
n_S <- 3; n_states_total <- 6

# 2. Format Counts (r) and Totals (n)
r <- array(NA, dim = c(n_intervs, n_states_total, n_states_total))
n_icd_total <- numeric(n_intervs)
n_noicd_total <- numeric(n_intervs)
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

for (i in seq_len(n_intervs)) {
  # Counts matrix
  r[i, , ] <- rbind(
    cbind(transdat[[i]]$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, transdat[[i]]$low_risk[, 1:n_S])
  )
  # Row Totals (Required for dmulti)
  n_icd_total[i]   <- sum(transdat[[i]]$ICD[1, 1:n_S])
  n_noicd_total[i] <- sum(transdat[[i]]$low_risk[1, 1:n_S])
}

# 3. Define Hyperpriors
# ---------------------
# Mean Log-Odds: -3 corresponds to prob ~4.7% (exp(-3)/(1+exp(-3)))
prior_mu <- c(-3, -3)

# Precision: Weakly informative
prior_tau <- diag(0.01, 2)

# Wishart (Covariance) Prior:
# Identity matrix scaled, with df=2 (minimally informative)
prior_R <- diag(0.1, 2)
prior_k <- 2

dataJags <- list(
  r = r,
  n_int = n_intervs,
  n_s = n_states_total,
  n_icd_total = n_icd_total,
  n_noicd_total = n_noicd_total,
  from_shock = c(1,0,0,0,0,0),
  from_ICD_death = c(0,0,1,0,0,0),
  prior_mu = prior_mu,
  prior_tau = prior_tau,
  prior_R = prior_R,
  prior_k = prior_k
)

filein <- "BUGS/model_bivariate_logit.txt"
params <- c("lambda", "rho_icd", "rho_noicd") # Monitor correlations!

# 4. Inits
inits <- function() {
  list(
    # Random initial values for the logits (around -3)
    beta_icd   = matrix(rnorm(n_intervs * 2, -3, 0.5), nrow = n_intervs, ncol = 2),
    beta_noicd = matrix(rnorm(n_intervs * 2, -3, 0.5), nrow = n_intervs, ncol = 2),

    # Init Hyperpriors
    mu_icd = c(-3, -3),
    mu_noicd = c(-3, -3),
    Tau_icd = diag(1, 2),
    Tau_noicd = diag(1, 2)
  )
}

# 5. Run Model (DIC = TRUE is okay here as it's a hierarchical fit)
tryCatch({
  mm1 <- jags(
    data = dataJags,
    inits = inits,
    parameters.to.save = params,
    model.file = filein,
    n.chains = 2,
    n.iter = 20000,
    n.burnin = 5000,
    n.thin = 5,
    DIC = TRUE
  )

  save(mm1, file = here::here("data/jags_output_all_strat.RData"))
  print("Model finished successfully.")
  print(mm1) # Check rho_icd and rho_noicd!

}, error = function(e) message(e))
