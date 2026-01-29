
# CE Markov model script
# main_ce_analysis_tunnel.R
# N Green, UCL
#
# using BUGS transition probabilities
# compare ICD vs non-ICD for different
# stratifications determined by
# 5 year risk prediction in Cox model


library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(ggplot2)
# library(HCM.SCD.CEanalysis)


########
# data #
########

# posterior samples of transition probabilities
# lambda.0: observed
# lambda.1: risk > 6%
# lambda.2: risk > 4%
# num [1:1000, 1:6, 1:6] = [samples, state, state]
data("jags_output")
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

# individual patient data
data("ipd_risk")

# risk score model fit
source(here::here("R/omahony_data.R"))

# model states:
# 1. HCM with ICD
# 2. shock
# 3. all-cause mortality
# 4. HCM (no ICD)
# 5. SCD
# 6. all-cause mortality

# structure:
# (1.icd  2.shock 3.mortality)
# (4.~icd 5.scd   6.mortality)


# multiple scenarios
inp_tab <-
  readr::read_csv(
    here::here("data/scenario_input_values.csv"),
    col_names = TRUE, show_col_types = FALSE)
##test values:
# inp_tab <-
#   readr::read_csv(
#     here::here("tests/testthat/scenario_input_values.csv"), col_names = TRUE)


#####################
# simulation set-up

n_tunnel <- inp_tab$t_repl[1]
n_total <- nrow(ipd_risk)
risk_thresh <- c(0.04, 0.06)
interv_names <- c("risk4", "risk6", "obs")
n_rules <- length(interv_names)
n_sim <- n.sims
S <- 2*n_tunnel + 4  # number of states
J <- 12              # max time (year)


#####################
# probabilities

# convert from no tunnel state to matrix with tunnel states
# list of decision rules
# e.g. num [1:6, 1:6, 1, 1:1000]  (no tunnel states)
# with time-homogeneous dimension length one
probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
alpha <- 0
probs <- list(obs = probs_empty,
              risk6 = probs_empty,
              risk4 = probs_empty)

for (i in 1:n_sim) {
  probs$obs[,,1,i] <-
    lambda_tunnel(lambda.0[i,,],
                  n_tunnel = n_tunnel,
                  alpha)
  probs$risk6[,,1,i] <-
    lambda_tunnel(lambda.1[i,,],
                  n_tunnel = n_tunnel,
                  alpha)
  probs$risk4[,,1,i] <-
    lambda_tunnel(lambda.2[i,,],
                  n_tunnel = n_tunnel,
                  alpha)
}

save(probs, file = "data/probs.RData")


######################################################
# sample init ICD pop using risk score distribution

# assume Normal distribution
# with frequentist mean and sd
# to approximate posterior

n_init_fixed <- c(risk4 = 2569, risk6 = 3134, obs = 3113)

norm_params <- tranformLN(hr)
coeff_samples <- rcoeff(norm_params, n_sim)

rf_names <- c("age", "mwt", "mwt2", "la", "maxlvotg",
              "fhxscd", "nsvt", "syncope")
rf <- ipd_risk[, rf_names]
r_init <- list()

for (j in risk_thresh) {

  p_risk <-
    pop_sample_by_age(j, rf, tmax = 2, coeff_samples)

  # move time to upper level
  # sample to lower level
  # extract only start time pop
  r_init[[as.character(j)]] <-
    map_dbl(purrr::transpose(p_risk)[[2]], length)
}

names(r_init) <- interv_names[1:2]

# append fixed obs value
##TODO: should/could we include some uncertainty?
r_init[["obs"]] <- n_init_fixed["obs"]

## suspend random init pop
## check: set as fixed point values
r_init[["risk4"]] <- n_init_fixed["risk4"]
r_init[["risk6"]] <- n_init_fixed["risk6"]


#################################################
# start state populations model input
# array with only two non-zero starting states

n_init <- list()
pop <- list()

for (i in interv_names) {
  n_init[[i]] <- matrix(0, nrow = length(r_init[[i]]), ncol = S)
  n_init[[i]][, 1] <- n_total - r_init[[i]]
  n_init[[i]][, S - 2] <- r_init[[i]]

  # include empty time dimension
  pop[[i]] <- array(NA, dim = c(S, J, n_sim))
  pop[[i]][, 1, ] <- t(n_init[[i]])
}

save(pop, file = "data/pop.RData")

# # summary statistics for paper
# quantile(n_init[[1]], probs = c(0.025, 0.5, 0.975))
# quantile(r_init[[2]]/3672, probs = c(0.025, 0.5, 0.975))
# quantile((3672 - r_init[[1]])/3672, probs = c(0.025, 0.5, 0.975))


#############
# run model #
#############

res <- list()

for (i in 1:nrow(inp_tab)) {

  inp <- inp_tab[i, ]

  res[[i]] <-
    ce_sim(pop,
           probs,
           inp,
           rc_unit,
           re_unit,
           rc_init,
           re_init,
           distn = FALSE,
           labels = interv_names)
}

save(res, file = here::here("data/res.RData"))


##########
# output #
##########

x11()
par(mfrow = c(3,4))

# for (i in 1) {
  for (i in 2:nrow(inp_tab)) {
  he <-
    bcea(-res[[i]]$e_total/n_total, -res[[i]]$c_total/n_total,
         ref = 3,
         interventions = interv_names)

  ceplane.plot(he, title = letters[i])
  # ceac.plot(he, title = letters[i])

  ## check
  # xx <- colMeans(res[[i]]$e_total)/n_total
  # e_means <- xx - xx["obs"]
  # yy <- colMeans(res[[i]]$c_total)/n_total
  # c_means <- yy - yy["obs"]
  # plot(e_means, c_means, xlim = c(-0.4,0.4), ylim = c(-2000, 6000))
  # text(e_means, c_means, labels = i)
}


# table
summary(he)
BCEA::tabulate_means(he)
apply(he$e, 2, mean)
apply(he$c, 2, mean)

save(he, file = here::here("data/bcea_data.RData"))


# proportions shocked

sum(res$pop$risk6[2, , 1])/filter(init, risk == "risk_over_6", status == "ICD")$n
# 0.229149
sum(res$pop$obs[2, , 1])/filter(init, risk == "icd", status == "ICD")$n
# 0.1866873
sum(res$pop$risk4[2, , 1])/filter(init, risk == "risk_over_4", status == "ICD")$n


########################################
# initial pop distribution histograms

for (i in 1:n_rules) {
  # png(file = glue::glue("images/init_hist_{risk_threshold[i]}.png"))
  hist(r_init[[i]], breaks = 25,
       xlim = c(500, 4000),
       main = "",
       xlab = "HCM starting population",
       freq = FALSE,
       border = "lightgrey")
  abline(v = mean(r_init[[i]]), col = "blue", lwd = 2)
  abline(v = n_init_fixed[i], col = "red", lwd = 2)

  d <- density(r_init[[i]], adjust = 2, to = 3672)
  lines(d, col = "red")
  # dev.off()
}


