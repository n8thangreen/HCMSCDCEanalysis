# State Population Counts Plot Over Time (Tunnel-State Model)
# Generates stacked percentage bar plots for TV-ICD and S-ICD

library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(purrr)
library(readr)
library(here)

devtools::load_all(here::here())

# 1. Load data
load(here::here("data/transdat.RData"))
data("ipd_risk")

load(here::here("data/jags_output.RData"))

lambda0_local <- mm1$BUGSoutput$sims.list$lambda.0
lambda1_local <- mm1$BUGSoutput$sims.list$lambda.1
lambda2_local <- mm1$BUGSoutput$sims.list$lambda.2
lambda3_local <- mm1$BUGSoutput$sims.list$lambda.3

n_sim <- mm1$BUGSoutput$n.sims
n_total <- nrow(ipd_risk)
J <- 30 # Time horizon (years)

inp_tab <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  col_names = TRUE, show_col_types = FALSE
)

init_counts <- get_init_counts(ipd_risk)
lambda_list <- list(obs = lambda0_local, risk6 = lambda1_local, risk4 = lambda2_local, no_icd = lambda3_local)

run_sim_and_get_pop <- function(device) {
  lbl <- tolower(device)
  inp <- inp_tab %>%
    dplyr::filter(Scenario_Label == lbl) %>%
    dplyr::slice(1)

  interv_names <- c("risk4", "risk6", "obs", "no_icd")
  n_tunnel <- inp$t_repl_max

  pop <- create_pop_array(init_counts, J, n_sim, n_tunnel = n_tunnel, interv_names = interv_names)
  p_repl <- sample_repl_probs(n_sim, inp$repl_median, inp$repl_shape, n_tunnel)
  probs <- create_probs_array(lambda_list, n_tunnel = n_tunnel, alpha = 0, p_repl = p_repl)

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

  # Extract and aggregate tunnel state populations into the 6 clinical health states
  plot_dat <- list()
  for (i in interv_names) {
    pop_arr <- out$pop[[i]] # (S, J, n_sim)
    mean_pop <- apply(pop_arr, c(1, 2), mean) # (S, J)

    agg_pop <- matrix(0, nrow = 6, ncol = J)
    agg_pop[1, ] <- colSums(mean_pop[1:n_tunnel, , drop = FALSE])
    agg_pop[2, ] <- colSums(mean_pop[(n_tunnel + 1):(2 * n_tunnel), , drop = FALSE])
    agg_pop[3, ] <- mean_pop[2 * n_tunnel + 1, ]
    agg_pop[4, ] <- mean_pop[2 * n_tunnel + 2, ]
    agg_pop[5, ] <- mean_pop[2 * n_tunnel + 3, ]
    agg_pop[6, ] <- mean_pop[2 * n_tunnel + 4, ]

    mean_df <- as.data.frame(agg_pop)
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

plot_and_save <- function(total_dat, device) {
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
      title = paste("State Population Over Time -", device, "(Tunnel Model)"),
      fill = "Health State"
    )

  filename <- paste0("images/state_pop_over_time_", tolower(device), "_tunnel.png")
  ggsave(facet_plot, file = filename, dpi = 320, width = 12, height = 6)
  cat("Saved plot to:", filename, "\n")
  return(facet_plot)
}

for (dev in c("TVICD", "SICD")) {
  cat("\nRunning simulation for:", dev, "\n")
  total_dat <- run_sim_and_get_pop(dev)
  p <- plot_and_save(total_dat, dev)
  if (dev == "TVICD") {
    ggsave(p, file = "images/state_pop_over_time.png", dpi = 320, width = 12, height = 6)
  }
}
