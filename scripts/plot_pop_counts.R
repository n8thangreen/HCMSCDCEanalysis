# state pop counts plot over time
# stacked percentage bar plot

library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(purrr)
library(readr)
library(here)

# Load package
devtools::load_all(here::here())

# 1. Load data
load(here::here("data/transdat.RData"))
data("ipd_risk")

# Load BUGS output directly (adhering to package rules)
load(here::here("data/jags_output.RData"))

lambda0_local <- mm1$BUGSoutput$sims.list$lambda.0
lambda1_local <- mm1$BUGSoutput$sims.list$lambda.1
lambda2_local <- mm1$BUGSoutput$sims.list$lambda.2
lambda3_local <- mm1$BUGSoutput$sims.list$lambda.3

n_sim <- mm1$BUGSoutput$n.sims
n_total <- nrow(ipd_risk)
J <- 30 # Time horizon for base-case analysis

# Read scenario parameter inputs
inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

# Setup initial counts
init <- inits_by_interv(ipd_risk, group_names = c("icd", "risk_over_6", "risk_over_4"))
init_counts <- list(
  obs = list(
    ICD    = dplyr::filter(init, risk == "icd", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "icd", status == "low_risk")$n
  ),
  risk6 = list(
    ICD    = dplyr::filter(init, risk == "risk_over_6", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "risk_over_6", status == "low_risk")$n
  ),
  risk4 = list(
    ICD    = dplyr::filter(init, risk == "risk_over_4", status == "ICD")$n,
    no_ICD = dplyr::filter(init, risk == "risk_over_4", status == "low_risk")$n
  ),
  no_icd = list(
    ICD    = 0,
    no_ICD = n_total
  )
)

# Function to run simulation and return aggregated populations
run_sim_and_get_pop <- function(device, model_type) {
  lbl <- tolower(device)
  inp <- inp_tab %>%
    dplyr::filter(Scenario_Label == lbl) %>%
    dplyr::slice(1)

  interv_names <- c("risk4", "risk6", "obs", "no_icd")

  if (model_type == "nontunnel") {
    S <- 6
    inp$t_repl <- 1 # Enforce t_repl = 1 for standard 6-state model
    n_tunnel <- 1

    n_init <- list(
      obs    = c(init_counts$obs$ICD,    0, 0, init_counts$obs$no_ICD,    0, 0),
      risk6  = c(init_counts$risk6$ICD,  0, 0, init_counts$risk6$no_ICD,  0, 0),
      risk4  = c(init_counts$risk4$ICD,  0, 0, init_counts$risk4$no_ICD,  0, 0),
      no_icd = c(init_counts$no_icd$ICD, 0, 0, init_counts$no_icd$no_ICD, 0, 0)
    )

    probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
    probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

    probs$obs[,, 1, ]    <- aperm(lambda0_local, c(2, 3, 1))
    probs$risk6[,, 1, ]  <- aperm(lambda1_local, c(2, 3, 1))
    probs$risk4[,, 1, ]  <- aperm(lambda2_local, c(2, 3, 1))
    probs$no_icd[,, 1, ] <- aperm(lambda3_local, c(2, 3, 1))

    pop <- list()
    for (i in interv_names) {
      pop[[i]] <- array(NA, dim = c(S, J, n_sim))
      pop[[i]][, 1, ] <- t(n_init[[i]])
    }
  } else {
    # tunnel state model
    n_tunnel <- inp$t_repl
    S <- 2 * n_tunnel + 4

    probs_empty <- array(NA, dim = c(S, S, 1, n_sim))
    probs <- list(obs = probs_empty, risk6 = probs_empty, risk4 = probs_empty, no_icd = probs_empty)

    alpha <- 0
    for (i in 1:n_sim) {
      probs$obs[,, 1, i]    <- lambda_tunnel(lambda0_local[i,,], n_tunnel = n_tunnel, alpha)
      probs$risk6[,, 1, i]  <- lambda_tunnel(lambda1_local[i,,], n_tunnel = n_tunnel, alpha)
      probs$risk4[,, 1, i]  <- lambda_tunnel(lambda2_local[i,,], n_tunnel = n_tunnel, alpha)
      probs$no_icd[,, 1, i] <- lambda_tunnel(lambda3_local[i,,], n_tunnel = n_tunnel, alpha)
    }

    pop <- list()
    for (i in interv_names) {
      pop[[i]] <- array(NA, dim = c(S, J, n_sim))
      pop[[i]][, 1, ] <- 0
      pop[[i]][1, 1, ] <- init_counts[[i]]$ICD
      pop[[i]][2 * n_tunnel + 2, 1, ] <- init_counts[[i]]$no_ICD
    }
  }

  out <- ce_sim(
    pop = pop,
    probs = probs,
    inp = inp,
    c_unit = rc_unit,
    e_unit = re_unit,
    c_init = rc_init,
    e_init = re_init,
    distn = FALSE,
    interv_names = interv_names
  )

  # Extract and aggregate state populations
  plot_dat <- list()
  for (i in interv_names) {
    pop_arr <- out$pop[[i]] # (S, J, n_sim)
    
    # Average across posterior simulation
    mean_pop <- apply(pop_arr, c(1, 2), mean) # (S, J)
    
    # Aggregate tunnel states to the 6 standard states if tunnel model
    if (model_type == "tunnel") {
      agg_pop <- matrix(0, nrow = 6, ncol = J)
      agg_pop[1, ] <- colSums(mean_pop[1:n_tunnel, , drop = FALSE])
      agg_pop[2, ] <- colSums(mean_pop[(n_tunnel + 1):(2 * n_tunnel), , drop = FALSE])
      agg_pop[3, ] <- mean_pop[2 * n_tunnel + 1, ]
      agg_pop[4, ] <- mean_pop[2 * n_tunnel + 2, ]
      agg_pop[5, ] <- mean_pop[2 * n_tunnel + 3, ]
      agg_pop[6, ] <- mean_pop[2 * n_tunnel + 4, ]
      mean_pop <- agg_pop
    }
    
    mean_df <- as.data.frame(mean_pop)
    state_names <- c("ICD_HCM", "Shock", "All_cause_ICD", "No_ICD_HCM", "SCD", "All_cause")
    mean_df$state <- state_names

    plot_dat[[i]] <- melt(mean_df, id.vars = "state", value.name = "pop", variable.name = "time") %>%
      group_by(state, time) %>%
      summarise(pop = sum(pop), .groups = "drop")

    plot_dat[[i]]$time <- as.numeric(gsub("V", "", plot_dat[[i]]$time))
    plot_dat[[i]]$state <- factor(plot_dat[[i]]$state,
                                  levels = c("SCD", "Shock", "All_cause_ICD", "All_cause", "ICD_HCM", "No_ICD_HCM"))
  }

  total_dat <- bind_rows(plot_dat, .id = "interv")
  total_dat$interv <- factor(total_dat$interv,
                              levels = c("no_icd", "obs", "risk6", "risk4"),
                              labels = c("No ICD Strategy", "Observed Strategy", "Cox 5-year SCD risk >6%", "Cox 5-year SCD risk >4%"))
  return(total_dat)
}

# Plot and save function
plot_and_save <- function(total_dat, device, model_type) {
  # Generate stacked percentage bar plot
  facet_plot <-
    ggplot(total_dat, aes(fill = state, y = pop, x = time)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_x_continuous(breaks = seq(0, J, by = 5)) +
    facet_grid(. ~ interv) +
    ylab("Population Proportion") +
    xlab("Time (years)") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    ) +
    labs(
      title = paste("State Population Over Time -", device, "(", model_type, "run)"),
      fill = "Health State"
    )

  filename <- paste0("images/state_pop_over_time_", tolower(device), "_", model_type, ".png")
  ggsave(facet_plot, file = filename, dpi = 320, width = 12, height = 6)
  cat("Saved plot to:", filename, "\n")
  return(facet_plot)
}

# Run all 4 combinations
combinations <- list(
  list(device = "TVICD", model_type = "nontunnel"),
  list(device = "TVICD", model_type = "tunnel"),
  list(device = "SICD", model_type = "nontunnel"),
  list(device = "SICD", model_type = "tunnel")
)

plots <- list()
for (comb in combinations) {
  cat("\nRunning simulation for:", comb$device, comb$model_type, "\n")
  total_dat <- run_sim_and_get_pop(comb$device, comb$model_type)
  p <- plot_and_save(total_dat, comb$device, comb$model_type)
  
  # Save default plot (TV-ICD nontunnel) copy as well
  if (comb$device == "TVICD" && comb$model_type == "nontunnel") {
    ggsave(p, file = "images/state_pop_over_time.png", dpi = 320, width = 12, height = 6)
    cat("Saved default plot copy to images/state_pop_over_time.png\n")
  }
}

# Generate the table for the default (TVICD nontunnel) run
cat("\nGenerating final populations table for default TVICD nontunnel...\n")
default_dat <- run_sim_and_get_pop("TVICD", "nontunnel")
plot_dat_list <- split(default_dat, default_dat$interv)

tab <-
  purrr::map(plot_dat_list,
             ~. |> group_by(state) |>
               summarise(pop = sum(pop))) |>
  plyr::join_all(by = "state") |>
  `colnames<-`(c("State", names(plot_dat_list)))

print(tab)
xtable::xtable(tab)
