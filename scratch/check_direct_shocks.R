library(dplyr)
devtools::load_all(here::here())

load(here::here("data/transdat.RData"))
load(here::here("data/jags_output_all_strat_direct.RData"))
data("ipd_risk")

n_total <- nrow(ipd_risk)
lambda_mean <- apply(mm1$BUGSoutput$sims.list$lambda, c(2, 3, 4), mean)

# copy background mortality
lambda_mean[, 2, 3] <- lambda_mean[, 1, 3]
lambda_mean[, 2, 1] <- 1 - lambda_mean[, 1, 3]

init <- inits_by_interv(ipd_risk, group_names = names(transdat))

for (nm in c("risk_over_4", "risk_over_6", "icd")) {
  # State 1 to 6
  pop0 <- c(
    filter(init, risk == nm, status == "ICD")$n,
    0, 0,
    filter(init, risk == nm, status == "low_risk")$n,
    0, 0
  )
  pop <- matrix(0, nrow=6, ncol=12)
  pop[, 1] <- pop0
  
  trans <- matrix(0, nrow=6, ncol=6)
  j_idx <- which(names(transdat) == nm)
  
  # transition matrix
  trans[1:3, 1:3] <- lambda_mean[j_idx, 1:3, 1:3]
  trans[4:6, 4:6] <- lambda_mean[j_idx, 4:6, 4:6]
  
  for (t in 2:12) {
    pop[, t] <- pop[, t-1] %*% trans
  }
  
  cat(nm, "shocks:", sum(pop[2, ]), "\n")
}
