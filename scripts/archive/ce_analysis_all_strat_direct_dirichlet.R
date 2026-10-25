# CE Markov model script
# ce_analysis_all_strat_direct_dirichlet.R
# N Green, UCL
#
# using BUGS transition probabilities
# compare ICD vs non-ICD for all
# stratifications determined by
# 5-year risk prediction in Cox model
# and number of risk factors

library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(readr)
library(here)
devtools::load_all(here::here())

########
# data #
########

# posterior samples of transition probabilities
# WinBUGS output

load(here::here("data/jags_output_all_strat_direct.RData"))

# Extract the lambda array: [Simulations, Interventions, From_State, To_State]
lambda_samples <- mm1$BUGSoutput$sims.list$lambda

load(here::here("data/transdat.RData"))

# individual patient data
# original extract
data("ipd_risk")

interv_names <- names(transdat)
n_intervs <- length(interv_names)
n_sim <- dim(lambda_samples)[1]
J <- 30

inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

init_counts <- inits_by_interv(ipd_risk, group_names = interv_names) |>
  dplyr::mutate(state = "hcm")

init_full <- expand.grid(
  risk = interv_names,
  status = c("ICD", "low_risk"),
  state = c("hcm", "scd/shock", "death")) |>
  dplyr::arrange(risk, status) |>
  dplyr::left_join(init_counts, by = c("risk", "status", "state")) |>
  tidyr::replace_na(list(n = 0))

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


#################
# probabilities #
#################

state_names <- c("hcm", "scd/shock", "death")
S <- length(state_names)
interv_names <- names(transdat)
# Note: The model uses 2*S = 6 states total (3 for ICD arm, 3 for No-ICD arm)

# add time homogeneous dimension
probs <- list()

for (j in seq_len(n_intervs)) {
  # 2. Format for CE Simulation
  # ce_sim expects: [From, To, Time, Sim]
  # We extract the slice for intervention 'j': lambda_samples[, j, , ]
  # This slice is [Sim, From, To] -> we aperm to [From, To, Sim]

  probs[[interv_names[j]]] <- array(NA, dim = c(2*S, 2*S, 1, n_sim))

  # Permute dimensions: [Sim, From, To] -> [From, To, 1, Sim]
  # (Sim is dim 1, From is dim 2, To is dim 3)
  probs[[interv_names[j]]][,,1,] <- aperm(lambda_samples[, j, , ], c(2, 3, 1))
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
  inp$t_repl_max <- 1     # no tunnel states

  res[[i]] <-
    ce_sim(pop, probs, inp,
           rc_unit,
           re_unit,
           rc_init,
           re_init,
           distn = FALSE,
           interv_names)
}

save(res, file = here::here("data/res_dsa_all_strat_direct.RData"))


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
          line = list(color = 7:1, type = c("solid", "solid", rep("dashed", 4), "solid")) +
  geom_vline(xintercept = c(20000, 30000), color = "grey")


png("images/contour_30yr_all_strat_direct.png")
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

ggsave(filename = "images/ceac_30yr_all_strat_direct.png", plot = gg, dpi = 640)

## check
# delta_c <- res$c_total - res$c_total[, "obs"]
# delta_e <- res$e_total - res$e_total[, "obs"]

