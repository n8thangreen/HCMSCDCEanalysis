# Cost-Effectiveness Comparative Plotting Script (All-Strat)
#
# Load saved S-ICD vs TV-ICD comparative results and generate the forest and scatter plots.

library(ggplot2)
library(ggrepel)
library(here)

# 1. Load data
rdata_file <- here::here("data/ce_analysis_all_strat_s-icd_results.RData")
if (!file.exists(rdata_file)) {
  stop("BCEA all-strat results file 'data/ce_analysis_all_strat_s-icd_results.RData' not found. Please run scripts/ce_analysis_all_strat_s-icd.R first to generate the comparative simulation output.")
}

load(rdata_file)
cat("Loaded comparative all-strat results. WTP Threshold:", WTP_THRESHOLD, "\n")

# Ensure images folder exists
dir.create(here::here("images"), showWarnings = FALSE)

# 2. Generate and save Incremental Net Benefit (INB) Forest Plot
cat("Generating S-ICD vs TV-ICD INB Forest Plot...\n")
gg_forest <- ggplot(ce_data, aes(x = INMB_Mean, y = reorder(Description, INMB_Mean))) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = INMB_Lower, xmax = INMB_Upper), height = 0.2, color = "grey40") +
  geom_point(size = 4, color = "#2c3e50") +
  geom_text(aes(label = round(INMB_Mean, 0)), vjust = -1.2, size = 4, fontface = "bold") +
  labs(
    title = "Incremental Net Benefit: S-ICD vs. TV-ICD",
    subtitle = "Comparing across 8 Risk Stratification Scenarios",
    x = paste0("INMB (£) at WTP £", WTP_THRESHOLD),
    y = "Risk Strategy",
    caption = "Positive values favor S-ICD (Subcutaneous). Negative values favor TV-ICD (Transvenous). Error bars are 95% CrI."
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(face = "italic", size = 12, margin = margin(b = 15)),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = here::here("images/inmb_strategies_forest.png"),
  plot = gg_forest,
  width = 8,
  height = 6,
  dpi = 320
)

# 3. Generate and save Cost-Effectiveness Plane (Strategy Scatter)
cat("Generating S-ICD vs TV-ICD Cost-Effectiveness Plane plot...\n")
gg_plane <- ggplot(ce_data, aes(x = Delta_E_Mean, y = Delta_C_Mean, label = Description)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_abline(intercept = 0, slope = WTP_THRESHOLD, linetype = "dashed", color = "grey50") +
  geom_point(size = 4, color = "#e74c3c", alpha = 0.8) +
  geom_text_repel(size = 4, fontface = "bold") +
  labs(
    title = "Cost-Effectiveness Plane (S-ICD vs. TV-ICD)",
    subtitle = "Mean outcomes for each Risk Strategy",
    x = "Incremental QALYs (S-ICD - TV-ICD)",
    y = "Incremental Cost (£) (S-ICD - TV-ICD)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(face = "italic", size = 12, margin = margin(b = 15)),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = here::here("images/ce_plane_strategies.png"),
  plot = gg_plane,
  width = 8,
  height = 6,
  dpi = 320
)
