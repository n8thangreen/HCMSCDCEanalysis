
# CE Markov model script
# main_ce_analysis.R
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
library(readr)
# library(HCM.SCD.CEanalysis)


########
# data #
########

# posterior samples of transition probabilities
# WinBUGS output
#
# lambda.0: observed
# lambda.1: risk > 6%
# lambda.2:
# lambda.3:
#
# num [1:1000, 1:6, 1:6] = [samples, state, state]
data("jags_output")
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

# individual patient data
# original extract
data("ipd_risk")

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


##########
# set-up #
##########

# start state populations
# point values
init <-
  inits_by_interv(
    ipd_risk,
    group_names = c("icd", "risk_over_6", "risk_over_4"))

n_init <- list()

n_init$obs <- c(
  filter(init, risk == "icd", status == "ICD")$n, 0, 0,
  filter(init, risk == "icd", status == "low_risk")$n, 0, 0)

n_init$risk6 <- c(
  filter(init, risk == "risk_over_6", status == "ICD")$n, 0, 0,
  filter(init, risk == "risk_over_6", status == "low_risk")$n, 0, 0)

n_init$risk4 <- c(
  filter(init, risk == "risk_over_4", status == "ICD")$n, 0, 0,
  filter(init, risk == "risk_over_4", status == "low_risk")$n, 0, 0)

n_init$no_icd <- c(0, 0, 0, 3672, 0, 0)

n_sim <- mm1$BUGSoutput$n.sims
S <- length(n_init$obs)  # number of states
J <- 30                  # max time (year)
#J <- 12                  # max time (year)
risk_thresh <- c(0.04, 0.06)
interv_names <- c("risk4", "risk6", "obs", "no_icd")


###################
# state values

## cost & health units
# multiple scenarios
# for sensitivity analysis
inp_tab <-
  readr::read_csv(
    here::here("data/scenario_input_values.csv"),
    col_names = TRUE, show_col_types = FALSE)


###################
# probabilities

# add time homogeneous dimension
probs_empty <- array(NA, dim = c(S, S, 1, n_sim))

probs <- list(obs = probs_empty,
              risk6 = probs_empty,
              risk4 = probs_empty,
              no_icd = probs_empty)

# include all-cause mortality from shock
for (i in 1:n_sim) {
  lambda.0[i, 2, 3] <- lambda.0[i, 1, 3]
  lambda.1[i, 2, 3] <- lambda.1[i, 1, 3]
  lambda.2[i, 2, 3] <- lambda.2[i, 1, 3]
  lambda.3[i, 2, 3] <- lambda.3[i, 1, 3]

  lambda.0[i, 2, 1] <- 1 - lambda.0[i, 2, 3]
  lambda.1[i, 2, 1] <- 1 - lambda.1[i, 2, 3]
  lambda.2[i, 2, 1] <- 1 - lambda.2[i, 2, 3]
  lambda.3[i, 2, 1] <- 1 - lambda.3[i, 2, 3]
}

# rearrange into 4 dim list
# num [1:6, 1:6, 1, 1:10000]
probs$obs[,,1,] <- aperm(lambda.0, c(2,3,1))
probs$risk6[,,1,] <- aperm(lambda.1, c(2,3,1))
probs$risk4[,,1,] <- aperm(lambda.2, c(2,3,1))
probs$no_icd[,,1,] <- aperm(lambda.3, c(2,3,1))


#################################################
# start state populations model input
# array with only two non-zero starting states

pop <- list()

for (i in interv_names) {
  # include empty time dimension
  pop[[i]] <- array(NA, dim = c(S, J, n_sim))
  pop[[i]][, 1, ] <- t(n_init[[i]])
}


#############
# run model #
#############

res <- list()

# over scenarios
for (i in 1:1){ #nrow(inp_tab)) {
  print(i)
  inp <- inp_tab[i, ]
  inp$t_repl <- 1     # no tunnel states

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

save(res, file = here::here("data/res_dsa.RData"))


##########
# output #
##########

library(BCEA)
library(ggplot2)

# baseline scenario
he <-
  bcea(eff = -res[[1]]$e_total/3672,
       cost = -res[[1]]$c_total/3672,
       ref = 3,
       .comparison = c(1,2),
       interventions = c("Cox 5-year SCD risk >4%", "Cox 5-year SCD risk >6%", "Observed", "no_icd"))

ceplane.plot(he, currency = "£", graph = "ggplot", title = "")

contour2(he, currency = "£", graph = "ggplot", title = "", ylim = c(-700, 2000), pos = c(1,0)) +
  theme(text = element_text(size = 20))
ggsave(filename = "images/contour_main.png", width = 7, height = 7, dpi = 640)


ceac.plot(he, title = "", pos = c(1,1), graph = "ggplot", currency = "£") +
  theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0.5, col = "grey")
ggsave(filename = "images/ceac_main.png", width = 7, height = 7, dpi = 640)

png("images/contour_30yr.png")
contour2(he)
dev.off()


# side-by-side for paper
he_multi_base <-
  multi.ce(
    bcea(eff = res[[1]]$e_total/3672,
       cost = res[[1]]$c_total/3672,
       ref = 3,
       interventions = c("risk4", "risk6", "obs", "no_icd")))

he_multi_u_icd95 <-
  multi.ce(
    bcea(eff = res[[18]]$e_total/3672,
         cost = res[[18]]$c_total/3672,
         ref = 3,
         interventions = c("risk4", "risk6", "obs", "no_icd")))

png("images/ceac_simulataneous_u_icd_grid.png", width = 25, height = 12, units = "cm", res = 320)
par(mfrow = c(1,2))
ceac.plot(he_multi_base, pos = c(1,1), currency = "£", graph = "base", title = ""); title(sub = "(a) Utility in HCM ICD state (u_icd) of 0.9")
ceac.plot(he_multi_u_icd95, pos = c(1,1), currency = "£", graph = "base", title = ""); title(sub = "(b) Utility in HCM ICD state (u_icd) of 0.95")
dev.off()



# gg <- ceac.plot(he_multi, pos = c(1,1), currency = "£", graph = "ggplot")
# gg
# ggsave(filename = "images/ceac_30yr.png", plot = gg, dpi = 640)






#########
# checks
# delta_c <- res$c_total - res$c_total[, "obs"]
# delta_e <- res$e_total - res$e_total[, "obs"]

