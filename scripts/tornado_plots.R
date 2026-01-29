
# one-way DSA tornado plots


library(ggplot2)
library(ceplot)
library(purrr)
library(BCEA)
library(readr)
library(dplyr)
library(ceplot)


load(file = here::here("data/res_dsa.RData"))
scenario_input_values <- read_csv("inst/extdata/scenario_input_values.csv")

N <- 3672
k <- 25000

he <-
  map(res,
      ~bcea(-.$e_total/N,
            -.$c_total/N,
            ref = 3,
            interventions = c("risk4", "risk6", "obs", "no_icd")))

# # test
# ce_total <-
#   rbind(cbind(name = "cost", id = 1:18,
#               -map_dfr(map(res, "c_total"), colMeans)/N),
#         cbind(name = "health", id = 1:18,
#               k*map_dfr(map(res, "e_total"), colMeans)/N))
#
# ce_total |>
#   group_by(id) |>
#   summarise(nb_risk4 = sum(risk4),
#             nb_risk6 = sum(risk6),
#             nb_obs = sum(obs),
#             nb_no_icd = sum(no_icd)) |>
#   mutate(eib_risk4 = nb_risk4 - nb_obs,
#          eib_risk6 = nb_risk6 - nb_obs,
#          eib_no_icd = nb_no_icd - nb_obs)


###########
# ICER

# tornado_dat <- data.frame(inp_tab, ICER = map_dfr(he, "ICER"))
#
# # 6%
# s_analysis <-
#   model.frame(formula = ICER.risk6 ~ q_hcm + u_icd + u_shock + u_implant + c_icd,
#               data = tornado_dat)
# ggtorn_dat <-
#   create_tornado_data(
#     s_analysis,
#     baseline_input = c(q_hcm=0.88, u_icd=0.9, u_shock=0.875, u_implant=-0.048, c_icd=4666))
# ggplot_tornado(ggtorn_dat, baseline_output = 1863) + ylab("ICER")
#
# # 4%
# s_analysis <-
#   model.frame(formula = ICER.risk4 ~ q_hcm + u_icd + u_shock + u_implant + c_icd,
#               data = tornado_dat)
# ggtorn_dat <-
#   create_tornado_data(
#     s_analysis,
#     baseline_input = c(q_hcm=0.88, u_icd=0.9, u_shock=0.875, u_implant=-0.048, c_icd=4666))
# ggplot_tornado(ggtorn_dat, baseline_output = -70407) + ylab("ICER")


###############
# EIB

eib_basecase <- he[[1]]$eib[he[[1]]$k == 25000, ]
basecase_inputs <- c(q_hcm = 0.88, u_icd = 0.9, u_shock = 0.875,
                     u_implant = -0.048, c_icd = 4666)

tornado_dat <-
  data.frame(scenario_input_values,
             EIB = map_dfr(he, function(x) x$eib[x$k == 25000, ]))

# 6% risk threshold
s_analysis <-
  model.frame(formula = EIB.risk6 ~ q_hcm + u_icd + u_shock + u_implant + c_icd,
              data = tornado_dat)

ggtorn_dat <-
  ceplot:::create_tornado_data(
    s_analysis,
    baseline_input = basecase_inputs)

torn_6 <- ceplot:::ggplot_tornado.tornado(
  dat = ggtorn_dat,
  baseline_output = eib_basecase[["risk6"]],
  annotate_nudge = 110) + ylab(paste0("EIB (", enc2utf8("\u00A3"), ")")) +
  labs(caption = "(b) > 6% 5-year SCD risk threshold Cox algorithm") +
  theme(plot.caption = element_text(hjust = 0, size = 12))


# 4% risk threshold
s_analysis <-
  model.frame(formula = EIB.risk4 ~ q_hcm + u_icd + u_shock + u_implant + c_icd,
              data = tornado_dat)

ggtorn_dat <-
  ceplot:::create_tornado_data(
    s_analysis,
    baseline_input = basecase_inputs)

torn_4 <- ceplot:::ggplot_tornado.tornado(
  dat = ggtorn_dat,
  baseline_output = eib_basecase[["risk4"]],
  annotate_nudge = 1000) + ylab(paste0("EIB (", enc2utf8("\u00A3"), ")")) +
  labs(caption = "(a) > 4% 5-year SCD risk threshold Cox algorithm") +
  theme(plot.caption = element_text(hjust = 0, size = 12))


library(grid)
library(gridExtra)
torn_grid <- grid_arrange_shared_legend(torn_4, torn_6, position = "right")

# save
ggsave(torn_6, filename = "images/tornado_plot_6_EIB.png",
       width = 4, height = 4, dpi = 640, scale = 1.5)

ggsave(torn_4, filename = "images/tornado_plot_4_EIB.png",
       width = 4, height = 4, dpi = 640, scale = 1.5)

ggsave(torn_grid, filename = "images/tornado_plots_EIB.png",
       width = 11, height = 6, dpi = 640)

