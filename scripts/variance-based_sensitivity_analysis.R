
# variance-based sensitivity analysis
# random sampling ce input parameter values


library(glue)
library(purrr)
library(dplyr)
library(BCEA)

rpoint <- function(n, x) x #rep(x, n)

input_distns <- readr::read_csv("data/input_distns.csv")
# for testing with no parameter uncertainty
# input_distns <- readr::read_csv("data/input_distns_test.csv")
load("data/pop.RData")
load("data/probs.RData")

# reduce number of iterations
# determined by the size of pop
pop <- map(pop, function(x) x[, , 1:100, drop = FALSE])

# replace all posterior samples with means
# a way to remove transition probability randomness
mean_probs <-
  list(risk4 = apply(probs$risk4, c(1,2,3), mean) %>%
         replicate(n = dim(pop$risk4)[3],
                   expr = .,
                   simplify = "array"),
       risk6 = apply(probs$risk6, c(1,2,3), mean) %>%
         replicate(n = dim(pop$risk6)[3],
                   expr = .,
                   simplify = "array"),
       obs = apply(probs$obs, c(1,2,3), mean) %>%
         replicate(n = dim(pop$obs)[3],
                   expr = .,
                   simplify = "array"))

## pre-sample all input values
# n_samples <- 10
# sample_fns <- glue("r{input_distns$distn}({n_samples}, {input_distns$params})")
# samples <- map(sample_fns, ~eval(parse(text = .)))
# model_ce_inputs <- expand.grid(samples)


## or on-the-fly during run-time

num_xi <- 100  # sample size
PSA_param_names <- c("q_hcm", "u_icd", "u_shock", "c_compl",
                     "c_shock", "c_icd", "c_appt")
interv_names <- c("risk4", "risk6", "obs")
wtp <- 25000
n_total <- 3672
res <- NULL
out <- list()

for (i in PSA_param_names) {
  print(paste0(i, ":"))

  out[[i]] <- list()
  row_idx <- which(input_distns$variable == i)
  sample_fn <-
    glue("r{input_distns$distn[row_idx]}(1, {input_distns$params[row_idx]})")

  for (j in seq_len(num_xi)) {
    print(j)
    # fix single input
    inp <- input_distns
    inp$distn[row_idx] <- "point"
    inp$params[row_idx] <- eval(parse(text = sample_fn))

    out[[i]][[j]] <-
      ce_sim(pop,
             mean_probs,
             # probs,
             inp,
             rc_unit,
             re_unit,
             rc_init,
             re_init,
             distn = TRUE,
             labels = interv_names)
    he <-
      bcea(-out[[i]][[j]]$e_total/n_total,
           -out[[i]][[j]]$c_total/n_total,
           ref = 3,
           interventions = interv_names)

    # choose output statistic
    res <-
      rbind.data.frame(res, c(i, he$eib[he$k == wtp, ]))
  }
  names(res) <- c("param", "risk4", "risk6")
}

# total variance (single run)
out_tot <- ce_sim(pop,
                  # mean_probs,
                  probs,
                  input_distns,
                  rc_unit,
                  re_unit,
                  rc_init,
                  re_init,
                  distn = TRUE,
                  labels = interv_names)
he_tot <-
  bcea(-out_tot$e_total/n_total,
       -out_tot$c_total/n_total,
       ref = 3,
       interventions = interv_names)

res_tot <- apply(he_tot$ib[he$k == wtp, , ], 2, var)

save(res, file = "data/variance_sa_res.RData")


##########################################
# variance-based SA statistics
# single parameter partial variance
# first-order effect

library(ggplot2)
library(reshape2)

plot_dat <-
  melt(res, id.vars = "param",
       variable.name = "interv",
       value.name = "eib") %>%
  group_by(param, interv) %>%
  summarise(var_partial = var(eib)) %>%
  ungroup() %>%
  mutate(var_total = res_tot[interv],
         var_explain = var_partial/var_total) %>%
  # hack to normalise because random
  # proportion can't be > 1
  group_by(interv) %>%
  mutate(var_explain = var_explain/sum(var_explain),
         interv = ifelse(interv == "risk4",
                         "4% threshold Cox model",
                         "6% threshold Cox model"))

ggplot(plot_dat,
       aes(x = param, y = var_explain, col = param, fill = param)) +
  geom_bar(stat = "identity") +
  facet_wrap(~interv) +
  theme_bw() + ylim(0,1) +
  theme(legend.position = "none") +
  xlab("Model parameter") +
  ylab("Variance explained")

ggsave(filename = "images/var_explain_barplot.png",
       width = 7, height = 5, dpi = 640)

# check total variation sums to one
plot_dat %>%
  group_by(interv) %>%
  summarise(sum(var_explain))


#######################
# EVPPI

inp_evi <-
  list(mat = `colnames<-`(out_tot$input_tab,
                          input_distns$variable),
       parameters = input_distns$variable)

# evppi(he_tot, param_idx = c(1,2,3), input = inp_evi$mat, N = 100)

info.rank(he_tot, inp = inp_evi, wtp = 25000)
info.rank(he_tot, inp = inp_evi, wtp = 25000, rel = FALSE)

ceac.plot(he_tot)
ceplane.plot(he_tot)

