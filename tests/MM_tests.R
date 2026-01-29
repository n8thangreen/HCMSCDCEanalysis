
# tests.R
# N Green


library(dplyr)
library(purrr)
library(BCEA)
library(reshape2)
library(HCM.SCD.CEanalysis)


n_init <- c(1, 0)
n_sim <- 2
S <- 2
J <- 12

lambda <- list(array(c(0.95, 0.85, 0.05, 0.15,
                       0.95, 0.85, 0.05, 0.15),       # drug A
                     dim = c(S, S, n_sim)),
               array(c(0.975, 0.95, 0.025, 0.05,
                       0.975, 0.95, 0.025, 0.05),     # drug B
                     dim = c(S, S, n_sim)))
c_unit <- list(c(50, 150),
               c(100, 200))
e_unit <- list(c(0.75, 0.73),
               c(0.75, 0.74))

res <-
  init_pop(n_init, n_sim, J) %>%
  ce_sim(lambda,
         c_unit,
         e_unit)


# sum across all time points

c <- map_dfc(res$cost, rowSums) %>% as.matrix()
e <- map_dfc(res$eff, rowSums) %>% as.matrix()

labels <- c("FP", "SFC")

m <- bcea(e, c,
          ref = 2,
          interventions = labels,
          Kmax = 300)
