# Variance-based sensitivity analysis (Tunnel-State Model)
# Random sampling CE input parameter values

library(glue)
library(purrr)
library(dplyr)
library(BCEA)
library(here)

devtools::load_all(here::here())

rpoint <- function(n, x) x

input_distns <- readr::read_csv(here::here("inst/extdata/input_distns.csv"), show_col_types = FALSE) %>%
  dplyr::filter(label == "tvicd")

data("ipd_risk")
load(here::here("data/jags_output.RData"))

n_sim <- 100
J <- 30
interv_names <- c("risk4", "risk6", "obs")

init_counts <- get_init_counts(ipd_risk)
pop <- create_pop_array(init_counts, J, n_sim, interv_names = interv_names)

lambda_list <- list(
  risk4 = mm1$BUGSoutput$sims.list$lambda.2[1:n_sim, , ],
  risk6 = mm1$BUGSoutput$sims.list$lambda.1[1:n_sim, , ],
  obs   = mm1$BUGSoutput$sims.list$lambda.0[1:n_sim, , ]
)
get_p <- function(v) as.numeric(input_distns$params[input_distns$variable == v])
p_repl <- sample_repl_probs(n_sim, get_p("repl_median"), get_p("repl_shape"))  # deterministic Weibull
probs <- create_probs_array(lambda_list, alpha = 0, p_repl = p_repl)

# Replace all posterior samples with means to isolate parameter uncertainty
mean_probs <- list(
  risk4 = apply(probs$risk4, c(1, 2, 3), mean) %>%
    replicate(n = n_sim, expr = ., simplify = "array"),
  risk6 = apply(probs$risk6, c(1, 2, 3), mean) %>%
    replicate(n = n_sim, expr = ., simplify = "array"),
  obs = apply(probs$obs, c(1, 2, 3), mean) %>%
    replicate(n = n_sim, expr = ., simplify = "array")
)

num_xi <- 100
PSA_param_names <- c("q_hcm", "u_icd", "u_shock", "c_compl",
                     "c_shock", "c_icd", "c_appt",
                     "p_compl", "p_repl_compl", "p_inf_yr", "p_dis_yr")
wtp <- 25000
n_total <- 3672
res <- NULL
out <- list()

for (i in PSA_param_names) {
  print(paste0(i, ":"))

  out[[i]] <- list()
  row_idx <- which(input_distns$variable == i)
  sample_fn <- glue("r{input_distns$distn[row_idx]}(1, {input_distns$params[row_idx]})")

  for (j in seq_len(num_xi)) {
    inp <- input_distns
    inp$distn[row_idx] <- "point"
    inp$params[row_idx] <- eval(parse(text = sample_fn))

    out[[i]][[j]] <- ce_sim(
      pop,
      mean_probs,
      inp,
      rc_unit,
      re_unit,
      rc_init,
      re_init,
      distn = TRUE,
      labels = interv_names
    )
    he <- bcea(
      -out[[i]][[j]]$e_total / n_total,
      -out[[i]][[j]]$c_total / n_total,
      ref = 3,
      interventions = interv_names
    )

    res <- rbind.data.frame(res, c(i, he$eib[he$k == wtp, ]))
  }
  names(res) <- c("param", "risk4", "risk6")
}

out_tot <- ce_sim(
  pop,
  probs,
  input_distns,
  rc_unit,
  re_unit,
  rc_init,
  re_init,
  distn = TRUE,
  labels = interv_names
)
he_tot <- bcea(
  -out_tot$e_total / n_total,
  -out_tot$c_total / n_total,
  ref = 3,
  interventions = interv_names
)

res_tot <- apply(he_tot$ib[he_tot$k == wtp, , ], 2, var)

save(res, file = here::here("data/variance_sa_res.RData"))

library(ggplot2)
library(reshape2)

plot_dat <-
  melt(res, id.vars = "param",
       variable.name = "interv",
       value.name = "eib") %>%
  mutate(eib = as.numeric(eib)) %>%
  group_by(param, interv) %>%
  summarise(var_partial = var(eib), .groups = "drop") %>%
  mutate(var_total = res_tot[as.character(interv)],
         var_explain = var_partial / var_total) %>%
  group_by(interv) %>%
  mutate(var_explain = var_explain / sum(var_explain),
         interv = ifelse(interv == "risk4",
                         "4% threshold Cox model",
                         "6% threshold Cox model"))

ggplot(plot_dat,
       aes(x = param, y = var_explain, col = param, fill = param)) +
  geom_bar(stat = "identity") +
  facet_wrap(~interv) +
  theme_bw() + ylim(0, 1) +
  theme(legend.position = "none") +
  xlab("Model parameter") +
  ylab("Variance explained")

ggsave(filename = here::here("images/var_explain_barplot.png"),
       width = 7, height = 5, dpi = 640)
