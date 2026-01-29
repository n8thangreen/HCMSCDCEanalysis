
# run BUGS model script
# cf BCEA book chapter 5
#
# estimate transition probabilities
# from discrete time counts
# using Dirichlet conjugacy directly
#
# transition matrix:
# 1 ICD      x x x 0 0 0
# 2 shock    1 0 0 0 0 0
# 3 death    0 0 1 0 0 0
# 4 low risk 0 0 0 x x x
# 5 scd      0 0 0 0 1 0
# 6 death    0 0 0 0 0 1


library(R2jags)
library(here)

# list of transition counts matrix
# for each scenario
data("transdat")

n_intervs <- length(transdat)

n_S <- nrow(transdat$icd$ICD)  # Number of active states (usually 3: ICD, Shock, Death)
n_states_total <- 2 * n_S      # Total states (6)

# Format Data for JAGS ---

# combine low and high risk
# into single transition matrix
# pad with 0s
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

## single diagonal counts matrices
r <- array(NA, dim = c(n_intervs, n_states_total, n_states_total))

for (i in seq_len(n_intervs)) {
  # Construct the 6x6 count matrix from the 3x3 blocks
  r[i, , ] <- rbind(
    cbind(transdat[[i]]$ICD[, 1:n_S], empty_mat),      # Row 1-3 (ICD start)
    cbind(empty_mat, transdat[[i]]$low_risk[, 1:n_S])) # Row 4-6 (No-ICD start)
}

# Define Priors & Parameters ---

alpha_prior <- c(1, 1, 1)  # Uninformative prior (Jeffreys would be 0.5)

dataJags <- list(
  r = r,
  n_int = n_intervs,
  n_s = n_states_total,
  from_shock = c(1,0,0,0,0,0),     # Fixed transitions from Shock
  from_ICD_death = c(0,0,1,0,0,0), # Fixed transitions from Death
  alpha_prior = alpha_prior
)

filein <- "BUGS/BUGS_code_model_all_strat_direct_dirichlet.txt"

params <- "lambda"  # probabilities

# Define Initial Values ---

inits <- function() {
  list(
    # Initialize probabilities for the 3 active states
    # rdirichlet is not base R, so we simulate using gammas normalized
    p_icd = structure(
      t(apply(matrix(1, n_intervs, 3), 1, function(x) {
        g <- rgamma(3, 1, 1)
        g / sum(g)
      })),
      .Dim = c(n_intervs, 3)
    ),

    p_no_icd = structure(
      t(apply(matrix(1, n_intervs, 3), 1, function(x) {
        g <- rgamma(3, 1, 1)
        g / sum(g)
      })),
      .Dim = c(n_intervs, 3)
    )
  )
}

# Run Model ---

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin) / 500)

mm1 <-
  jags(data = dataJags,
       inits = inits,
       parameters.to.save = params,
       model.file = filein,
       n.chains = 2,
       n.iter,
       n.burnin,
       n.thin,
       DIC = FALSE)

R2WinBUGS::attach.bugs(mm1$BUGSoutput)

save(mm1, file = here::here("data/jags_output_all_strat_direct.RData"))
saveRDS(lambda, file = here::here("data/lambda_direct.Rds"))

# print(mm1, digits = 3, intervals = c(0.025, 0.975))

