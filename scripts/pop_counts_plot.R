
# state pop counts plot over time
# stacked percentage bar plot


library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)


load(file = "data/res_dsa.RData")

plot_dat <- list()

for (i in names(res[[1]]$pop)) {

  # average across posterior simulation
  plot_dat[[i]] <-
    apply(res[[1]]$pop[[i]], c(1,2), mean) %>%
    as.data.frame()

  state_names <-  c("ICD_HCM", "Shock", "All_cause_ICD", "No_ICD_HCM", "SCD", "All_cause")
  # state_names <-  c("ICD_HCM", "Shock", "All_cause", "No_ICD_HCM", "SCD", "All_cause")

  plot_dat[[i]]$state <- state_names

  plot_dat[[i]] <-
    melt(plot_dat[[i]],
         id.vars = "state",
         value.name = "pop",
         variable.name = "time") %>%
    group_by(state, time) %>%
    summarise(pop = sum(pop))

  # remove from time label
  plot_dat[[i]]$time <- as.numeric(gsub("V", "", plot_dat[[i]]$time))
  plot_dat[[i]]$state <-
    factor(plot_dat[[i]]$state,
           levels = c("SCD", "Shock", "All_cause_ICD", "All_cause", "ICD_HCM", "No_ICD_HCM"))
           # levels = c("SCD", "Shock", "All_cause", "ICD_HCM", "No_ICD_HCM"))

}

risk_scores <- as_labeller(c(
  'no_icd' = "No ICDs",
  'obs' = "Observed data",
  'risk4' = "Cox 5-year risk >4%",
  'risk6'="Cox 5-year risk >6%"))

total_dat <- bind_rows(plot_dat, .id = "interv")

facet_plot <-
  ggplot(total_dat, aes(fill = state, y = pop, x = time)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_x_continuous(breaks = 1:12) +
  facet_grid(. ~ interv, labeller = risk_scores) +
  ylab("Population") +
  xlab("Time (years)")

facet_plot

ggsave(facet_plot, file = "images/state_pop_over_time.png", dpi = 640, width = 10)


## table of final populations
tab <-
  purrr::map(plot_dat,
             ~. |> group_by(state) |>
               # filter(time == 12) |>
               summarise(pop = sum(pop))) |>
  plyr::join_all(by = "state") |>
  `colnames<-`(c("State", names(plot_dat)))
tab

xtable::xtable(tab)

