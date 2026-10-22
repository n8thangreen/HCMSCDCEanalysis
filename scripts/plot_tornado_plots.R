# One-way DSA Tornado Plots (TV-ICD and S-ICD Tunnel-State Models)

library(ggplot2)
library(ceplot)
library(purrr)
library(BCEA)
library(readr)
library(dplyr)
library(here)
library(grid)
library(gridExtra)

devtools::load_all(here::here())

load(file = here::here("data/res_dsa.RData"))
scenario_input_values <- readr::read_csv(
  here::here("inst/extdata/scenario_input_values.csv"),
  show_col_types = FALSE
)

N <- 3672
k <- 25000

col_names <- colnames(res[[1]]$e_total)
cols_to_compare <- if ("risk4" %in% col_names) {
  c("risk4", "risk6", "obs", "no_icd")
} else {
  c("risk_over_4", "risk_over_6", "icd", "no_icd")
}

he <- map(
  res,
  ~bcea(
    -.$e_total[, cols_to_compare] / N,
    -.$c_total[, cols_to_compare] / N,
    ref = 3,
    interventions = c("risk4", "risk6", "obs", "no_icd")
  )
)

tornado_dat_all <- data.frame(
  scenario_input_values,
  EIB = map_dfr(he, function(x) x$eib[x$k == k, ])
)

generate_device_tornado <- function(dev_label, suffix = "") {
  # Filter to the 13 core one-way DSA rows (baseline + 12 parameter bounds) for this device
  dev_rows <- tornado_dat_all %>%
    dplyr::filter(Scenario_Label == dev_label) %>%
    dplyr::slice(1:13)

  base_row <- dev_rows[1, ]
  basecase_inputs <- c(
    q_hcm     = base_row$q_hcm,
    u_icd     = base_row$u_icd,
    u_shock   = base_row$u_shock,
    u_implant = base_row$u_implant,
    c_icd     = base_row$c_icd,
    c_shock   = base_row$c_shock
  )

  # 6% risk threshold
  s_analysis_6 <- model.frame(
    formula = EIB.risk6 ~ q_hcm + u_icd + u_shock + u_implant + c_icd + c_shock,
    data = dev_rows
  )
  ggtorn_dat_6 <- ceplot:::create_tornado_data(
    s_analysis_6,
    baseline_input = basecase_inputs
  )
  torn_6 <- ceplot:::ggplot_tornado.tornado(
    dat = ggtorn_dat_6,
    baseline_output = base_row$EIB.risk6,
    annotate_nudge = 110
  ) +
    ylab(paste0("EIB (", enc2utf8("\u00A3"), ")")) +
    labs(caption = paste0("(b) > 6% 5-year SCD risk threshold (", toupper(dev_label), ")")) +
    theme(plot.caption = element_text(hjust = 0, size = 11))

  # 4% risk threshold
  s_analysis_4 <- model.frame(
    formula = EIB.risk4 ~ q_hcm + u_icd + u_shock + u_implant + c_icd + c_shock,
    data = dev_rows
  )
  ggtorn_dat_4 <- ceplot:::create_tornado_data(
    s_analysis_4,
    baseline_input = basecase_inputs
  )
  torn_4 <- ceplot:::ggplot_tornado.tornado(
    dat = ggtorn_dat_4,
    baseline_output = base_row$EIB.risk4,
    annotate_nudge = 500
  ) +
    ylab(paste0("EIB (", enc2utf8("\u00A3"), ")")) +
    labs(caption = paste0("(a) > 4% 5-year SCD risk threshold (", toupper(dev_label), ")")) +
    theme(plot.caption = element_text(hjust = 0, size = 11))

  torn_grid <- grid_arrange_shared_legend(torn_4, torn_6, position = "right")

  ggsave(torn_6, filename = paste0("images/tornado_plot_6_EIB", suffix, ".png"),
         width = 4, height = 4, dpi = 640, scale = 1.5)
  ggsave(torn_4, filename = paste0("images/tornado_plot_4_EIB", suffix, ".png"),
         width = 4, height = 4, dpi = 640, scale = 1.5)
  ggsave(torn_grid, filename = paste0("images/tornado_plots_EIB", suffix, ".png"),
         width = 11, height = 6, dpi = 640)
}

generate_device_tornado("tvicd", suffix = "")
generate_device_tornado("sicd", suffix = "_sicd")
