
# run BUGS model script
# cf BCEA book chapter 5
#
# estimate transition probabilities
# from discrete time counts
# using multinomial likelihood
#
# transition matrix:
# 1 ICD      x x x 0 0 0
# 2 shock    1 0 0 0 0 0
# 3 death    0 0 1 0 0 0
# 4 low risk 0 0 0 x x x
# 5 scd      0 0 0 0 1 0
# 6 death    0 0 0 0 0 1

# lambda.0: observed
# lambda.1: >6% risk group
# lambda.2: >4% risk group
# lambda.3: no icd


library(R2jags)

data("transdat")

data_obs <- transdat$icd
data_risk6 <- transdat$risk_over_6
data_risk4 <- transdat$risk_over_4
data_no_icd <- transdat$no_icd

n_S <- nrow(data_obs$ICD)  # number of states

# combine low and high risk
# into single transition matrix
# pad with 0s
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

## single diagonal counts matrices

# observed
r.0 <-
  rbind(
    cbind(data_obs$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, data_obs$low_risk[, 1:n_S]))

# risk >6
r.1 <-
  rbind(
    cbind(data_risk6$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, data_risk6$low_risk[, 1:n_S]))

# risk >4
r.2 <-
  rbind(
    cbind(data_risk4$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, data_risk4$low_risk[, 1:n_S]))

# no ICD
r.3 <-
  rbind(
    cbind(data_no_icd$ICD[, 1:n_S], empty_mat),
    cbind(empty_mat, data_no_icd$low_risk[, 1:n_S]))

# total from populations
n.0 <- unname(rowSums(r.0))
n.1 <- unname(rowSums(r.1))
n.2 <- unname(rowSums(r.2))
n.3 <- unname(rowSums(r.3))

scale <- 1                                 # level of informativeness for
alpha.0 <- c(rep(scale, n_S), rep(0, n_S)) # the Dirichlet prior
alpha.1 <- c(rep(0, n_S), rep(scale, n_S)) # with structural zeros

dataJags <-
  list(n.0 = n.0,
       n.1 = n.1,
       n.2 = n.2,
       n.3 = n.3,
       r.0 = r.0,
       r.1 = r.1,
       r.2 = r.2,
       r.3 = r.3,
       from_shock = c(1,0,0,0,0,0),
       from_ICD_death = c(0,0,1,0,0,0),
       alpha.0 = alpha.0,
       alpha.1 = alpha.1)

filein <- "BUGS/model_with_no_icd.txt"

# probabilities
params <- c("lambda.0",
            "lambda.1",
            "lambda.2",
            "lambda.3")

# intital values
inits <- function() {

  matgam.0 <- matrix(rgamma(2*2*n_S, scale, 1),
                     nrow = 2)
  # prohibited transitions
  matgam.0[1, (n_S+1):(2*n_S)] <- 0
  matgam.0[2, 1:n_S] <- 0
  sum.matgam.0 <- rowSums(matgam.0, na.rm = TRUE)

  matgam.1 <- matrix(rgamma(2*2*n_S, scale, 1),
                     nrow = 2)
  # prohibited transitions
  matgam.1[1, (n_S+1):(2*n_S)] <- 0
  matgam.1[2, 1:n_S] <- 0
  sum.matgam.1 <- rowSums(matgam.1, na.rm = TRUE)

  matgam.2 <- matrix(rgamma(2*2*n_S, scale, 1),
                     nrow = 2)
  # prohibited transitions
  matgam.2[1, (n_S+1):(2*n_S)] <- 0
  matgam.2[2, 1:n_S] <- 0
  sum.matgam.2 <- rowSums(matgam.2, na.rm = TRUE)

  matgam.3 <- matrix(rgamma(2*2*n_S, scale, 1),
                     nrow = 2)
  # prohibited transitions
  matgam.3[1, (n_S+1):(2*n_S)] <- 0
  matgam.3[2, 1:n_S] <- 0
  sum.matgam.3 <- rowSums(matgam.3, na.rm = TRUE)

  p.0 <- matgam.0/sum.matgam.0
  p.1 <- matgam.1/sum.matgam.1
  p.2 <- matgam.2/sum.matgam.2
  p.3 <- matgam.3/sum.matgam.3

  obs_mat <-
    rbind(p.0,
          matrix(NA,
                 ncol = 2*n_S,
                 nrow = 2*n_S - 2))

  risk6_mat <-
    rbind(p.1,
          matrix(NA,
                 ncol = 2*n_S,
                 nrow = 2*n_S - 2))

  risk4_mat <-
    rbind(p.2,
          matrix(NA,
                 ncol = 2*n_S,
                 nrow = 2*n_S - 2))

  no_icd_mat <-
    rbind(p.3,
          matrix(NA,
                 ncol = 2*n_S,
                 nrow = 2*n_S - 2))

  # rearrange rows
  new_order <- c(1,3,4,2,5,6)

  list(lambda.0 = obs_mat[new_order, ],
       lambda.1 = risk6_mat[new_order, ],
       lambda.2 = risk4_mat[new_order, ],
       lambda.3 = no_icd_mat[new_order, ])
}

n.iter <- 10000
n.burnin <- 5000
n.thin <- floor((n.iter - n.burnin)/500)

mm1 <-
  jags(data = dataJags,
       inits = inits,
       parameters.to.save = params,
       model.file = filein,
       n.chains = 2,
       n.iter,
       n.burnin,
       n.thin,
       DIC = TRUE)

R2WinBUGS::attach.bugs(mm1$BUGSoutput)

save(mm1, file = here::here("data/jags_output.RData"))
saveRDS(lambda.0, file = here::here("data/lambda0.Rds"))  # observed
saveRDS(lambda.1, file = here::here("data/lambda1.Rds"))  # > 6%
saveRDS(lambda.2, file = here::here("data/lambda2.Rds"))  # > 4%
saveRDS(lambda.3, file = here::here("data/lambda3.Rds"))  # no icd

# print(mm1, digits = 3, intervals = c(0.025, 0.975))

