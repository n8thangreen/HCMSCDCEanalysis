
# CE Markov model script
# all_strat_ce_analysis.R
# N Green, UCL
#
# using BUGS transition probabilities
# compare ICD vs non-ICD for all
# stratifications determined by
# 5 year risk prediction in Cox model
# and number of risk factors

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
# num [1:1000, 1:6, 1:6] = [samples, state, state]
load("data/jags_output_all_strat.RData")
R2WinBUGS::attach.bugs(mm1$BUGSoutput)

load("data/transdat.RData")

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

state_names <- c("hcm", "scd/shock", "death")
interv_names <- names(transdat)

# start state populations
init <-
  inits_by_interv(
    ipd_risk,
    group_names = interv_names) |>
  mutate(state = "hcm")

# fill in non-zero pops
init_full <-
  expand.grid(risk = interv_names,
              status = c("ICD", "low_risk"),
              state = state_names) |>
  arrange(risk, status) |>
  left_join(init) |>
  tidyr::replace_na(list(n = 0))

n_sim <- mm1$BUGSoutput$n.sims
S <- length(state_names)
n_intervs <- length(interv_names)
# J <- 30                  # max time (year)
J <- 12


###################
# state values

## cost & health units
# multiple scenarios
# for sensitivity analysis
inp_tab <-
  readr::read_csv(
    here::here("inst/extdata/scenario_input_values.csv"),
    col_names = TRUE, show_col_types = FALSE)


###################
# probabilities

# add time homogeneous dimension
probs <- list()

# include all-cause mortality from shock
# same as from hcm icd state
for (j in seq_len(n_intervs)) {

  for (i in seq_len(n_sim)) {
    lambda[i, j, 2, 3] <- lambda[i, j, 1, 3]
    lambda[i, j, 2, 1] <- 1 - lambda[i, j, 2, 3]
  }

  probs[[interv_names[j]]] <- array(NA, dim = c(2*S, 2*S, 1, n_sim))

  # rearrange into list of 4 dim arrays
  # num [1:6, 1:6, 1, 1:10000]
  probs[[interv_names[j]]][,,1,] <- aperm(lambda[, j, , ], c(2,3,1))
}


#################################################
# start state populations model input
# array with only two non-zero starting states

pop0 <- split(init_full$n, init_full$risk)
pop <- list()

for (i in interv_names) {
  # include empty time dimension
  pop[[i]] <- array(NA, dim = c(2*S, J, n_sim))
  pop[[i]][, 1, ] <- pop0[[i]]
}


#############
# run model #
#############

res <- list()

# over scenarios
# for (i in 1:1){
for (i in 1:nrow(inp_tab)) {
  print(i)
  inp <- inp_tab[i, ]
  inp$t_repl <- 1     # no tunnel states

  res[[i]] <-
    ce_sim(pop, probs, inp,
           rc_unit,
           re_unit,
           rc_init,
           re_init,
           distn = FALSE,
           interv_names)
}

save(res, file = here::here("data/res_dsa_all_strat.RData"))


##########
# output #
##########

library(BCEA)
library(ggplot2)

n_sample <- 3672

# 18: u_icd = 0.95
he <-
  bcea(eff = -res[[18]]$e_total/n_sample,
       cost = -res[[18]]$c_total/n_sample,
       ref = 7,
       interventions = interv_names)

ceplane.plot(he, graph = "ggplot2", point = list(color = 1:7), title = "")
ceac.plot(he, graph = "ggplot2", title = "", pos = c(1,1), currency = "£")#,
          line = list(color = 7:1, type = c("solid", "solid", rep("dashed", 4), "solid"))) +
  geom_vline(xintercept = c(20000, 30000), color = "grey")


png("images/contour_30yr_all_strat.png")
contour2(he)
dev.off()

he_multi <-
  multi.ce(
    bcea(eff = res[[18]]$e_total/n_sample,
       cost = res[[18]]$c_total/n_sample,
       ref = 7,
       interventions = interv_names))

gg <-
  ceac.plot(he_multi, graph = "ggplot2", title = "", pos = c(1,1), currency = "£",
            line = list(color = c(7:1,8), type = c("solid", "solid", rep("dashed", 4), "solid", "solid"))) +
  geom_vline(xintercept = c(20000, 30000), color = "grey")

gg

ggsave(filename = "images/ceac_30yr_all_strat.png", plot = gg, dpi = 640)

## check
# delta_c <- res$c_total - res$c_total[, "obs"]
# delta_e <- res$e_total - res$e_total[, "obs"]

