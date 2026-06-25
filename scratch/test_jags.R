# Test JAGS compilation
library(R2jags)
load(here::here("data/transdat.RData"))

data_obs <- transdat$icd
data_risk6 <- transdat$risk_over_6
data_risk4 <- transdat$risk_over_4
data_no_icd <- transdat$no_icd

n_S <- nrow(data_obs$ICD)
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

r.0 <- rbind(cbind(data_obs$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_obs$low_risk[, 1:n_S]))
r.1 <- rbind(cbind(data_risk6$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_risk6$low_risk[, 1:n_S]))
r.2 <- rbind(cbind(data_risk4$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_risk4$low_risk[, 1:n_S]))
r.3 <- rbind(cbind(data_no_icd$ICD[, 1:n_S], empty_mat),
             cbind(empty_mat, data_no_icd$low_risk[, 1:n_S]))

n.0 <- unname(rowSums(r.0))
n.1 <- unname(rowSums(r.1))
n.2 <- unname(rowSums(r.2))
n.3 <- unname(rowSums(r.3))

n_icd <- c(n.0[1], n.1[1], n.2[1], n.3[1])
n_noicd <- c(n.0[4], n.1[4], n.2[4], n.3[4])

r_icd_shock <- c(r.0[1, 2], r.1[1, 2], r.2[1, 2], r.3[1, 2])
r_icd_death <- c(r.0[1, 3], r.1[1, 3], r.2[1, 3], r.3[1, 3])
r_noicd_scd <- c(r.0[4, 5], r.1[4, 5], r.2[4, 5], r.3[4, 5])
r_noicd_death <- c(r.0[4, 6], r.1[4, 6], r.2[4, 6], r.3[4, 6])

mean_prior <- c(-5, -5)
prec_prior <- diag(2) * 0.1
R_wishart <- diag(2) * 0.1

print("n_icd:")
print(n_icd)
print("r_icd_shock:")
print(r_icd_shock)
print("r_icd_death:")
print(r_icd_death)

dataJags <- list(
  n_icd = n_icd,
  n_noicd = n_noicd,
  r_icd_shock = r_icd_shock,
  r_icd_death = r_icd_death,
  r_noicd_scd = r_noicd_scd,
  r_noicd_death = r_noicd_death,
  mean_prior = mean_prior,
  prec_prior = prec_prior,
  R_wishart = R_wishart,
  from_shock = c(1, 0, 0, 0, 0, 0),
  from_ICD_death = c(0, 0, 1, 0, 0, 0)
)

inits <- function() {
  list(
    p_icd_shock = rep(0.05, 4),
    p_noicd_scd = rep(0.01, 4),
    logit_p_death = matrix(rnorm(8, -5, 0.5), nrow = 4, ncol = 2),
    Omega_death = diag(2)
  )
}

mm1 <- jags(
  data = dataJags,
  inits = inits,
  parameters.to.save = c("lambda.0", "rate_icd_shock"),
  model.file = "BUGS/BUGS_code_best_model.txt",
  n.chains = 2,
  n.iter = 100,
  n.burnin = 10
)
print("Done compiles!")
