
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


library(R2jags)

data("transdat")

n_intervs <- length(transdat)

n_S <- nrow(transdat$icd$ICD)  # number of states

# combine low and high risk
# into single transition matrix
# pad with 0s
empty_mat <- matrix(0, nrow = n_S, ncol = n_S)

## single diagonal counts matrices
r <- array(NA, dim = c(n_intervs, 2*n_S, 2*n_S))
n <- NULL

for (i in seq_len(n_intervs)) {
  r[i, , ] <-
    rbind(
      cbind(transdat[[i]]$ICD[, 1:n_S], empty_mat),
      cbind(empty_mat, transdat[[i]]$low_risk[, 1:n_S]))

  # total from populations
  n <- rbind(n, unname(rowSums(r[i, , ])))
}

scale <- 1                                 # level of informativeness for
alpha.0 <- c(rep(scale, n_S), rep(0, n_S)) # the Dirichlet prior
alpha.1 <- c(rep(0, n_S), rep(scale, n_S)) # with structural zeros

dataJags <-
  list(n = n,
       r = r,
       n_int = n_intervs,
       n_s = 2*n_S,
       from_shock = c(1,0,0,0,0,0),
       from_ICD_death = c(0,0,1,0,0,0),
       alpha.0 = alpha.0,
       alpha.1 = alpha.1)

filein <- "BUGS/model_all_strat.txt"

# probabilities
params <- "lambda"

# initial transition probability values
inits <- function() {

  mat <- array(dim = c(n_intervs, 2*n_S, 2*n_S))

  for (i in seq_len(n_intervs)) {
    matgam <- matrix(rgamma(2*2*n_S, scale, 1),
                            nrow = 2)

    # prohibited transitions
    matgam[1, (n_S+1):(2*n_S)] <- 0
    matgam[2, 1:n_S] <- 0
    sum.matgam <- rowSums(matgam, na.rm = TRUE)

    p <- matgam/sum.matgam

    mat[i, , ] <-
      rbind(p, matrix(data = NA,
                      ncol = 2*n_S,
                      nrow = 2*n_S - 2))

    # rearrange rows
    new_order <- c(1,3,4,2,5,6)
    mat[i, , ] <- mat[i, new_order, ]
  }

  list("lambda" = mat)
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

save(mm1, file = here::here("data/jags_output_all_strat.RData"))
saveRDS(lambda, file = here::here("data/lambda.Rds"))

# print(mm1, digits = 3, intervals = c(0.025, 0.975))

