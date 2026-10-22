# Cost-Effectiveness Plotting Script
#
# Load saved cost-effectiveness results and generate standard decision plots.

library(BCEA)
library(ggplot2)
library(here)

# 1. Load data
rdata_file <- here::here("data/bcea_data.RData")

if (!file.exists(rdata_file)) {
  stop("BCEA data file 'data/bcea_data.RData' not found. Please run scripts/run_ce_analysis.R first to generate the simulation output.")
}

load(rdata_file)
cat("Loaded BCEA data for Device:", device, "| Model Structure:", model_type, "\n")

# Ensure images folder exists
dir.create(here::here("images"), showWarnings = FALSE)

# 2. Generate and save Cost-effectiveness Plane
cat("Generating Cost-effectiveness Plane plot...\n")
ce_plane_p <- ceplane.plot(he, currency = "£", graph = "ggplot", title = paste("CE Plane for", device))
ce_plane_p

ggsave(
  filename = here::here(paste0("images/ce_plane_", device, "_", model_type, ".png")),
  plot = ce_plane_p,
  width = 7,
  height = 7,
  dpi = 320
)

contour2(he)

# 3. Generate and save Cost-effectiveness Acceptability Curve (CEAC)
cat("Generating Cost-effectiveness Acceptability Curve (CEAC) plot...\n")
ceac_p <- ceac.plot(he, title = paste("CEAC for", device), graph = "ggplot", currency = "£") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5)) +
  geom_hline(yintercept = 0.5, col = "grey50", linetype = "dashed")
ceac_p

ggsave(
  filename = here::here(paste0("images/ceac_", device, "_", model_type, ".png")),
  plot = ceac_p,
  width = 7,
  height = 7,
  dpi = 320
)
