# test-ce_sim


## mock inputs

S <- 2      # number of states
n_sim <- 2
t_lim <- 10

# transition probabilities
# 1. always stay
# 2. always move
p <- array(c(1,0,
             0,1,
             0,1,
             1,0),
           dim = c(S, S, 1, n_sim))

# two interventions
probs = list(p, p)

# interv > state x time
c_unit <- function(S, t_lim) {
  function() {
    c_state_by_time <-
      matrix(1,
             ncol = t_lim,
             nrow = S,
             byrow = FALSE)

    list(c_state_by_time,
         c_state_by_time)
  }
}

# interv > state
c_init <- function() {
  function()
    list(c(1,1), c(1,1))}

e_unit <- function(S, t_lim) {
  function() {
    list(matrix(1,
                ncol = t_lim,
                nrow = S,
                byrow = FALSE),
         matrix(0,
                ncol = t_lim,
                nrow = S,
                byrow = FALSE))
  }
}

# interv > state x time x sim
pop <- list()
for (i in 1:2) {
  pop[[i]] <- array(NA, dim = c(S, t_lim, n_sim))
  pop[[i]][, 1, ] <- c(1,0)
}

test_that("basic inputs", {

  out <-
    ce_sim(pop = pop,
           probs = probs,
           c_unit = c_unit(S, t_lim),
           e_unit = e_unit(S, t_lim),
           c_init = c_init())

  ## state populations over time

  # stay in the starting state
  expect_equal(out$pop[[1]][, , 1],
               rbind(rep(1, t_lim),
                     rep(0, t_lim)))

  # move every time step
  expect_equal(out$pop[[1]][, , 2],
               rbind(rep(c(1,0), t_lim/2),
                     rep(c(0,1), t_lim/2)))

  ## state costs over time

  # same between simulations
  expect_equal(out$cost[[1]][1, ],
               out$cost[[1]][2, ])
  expect_equal(out$dcost[[1]][1, ],
               out$dcost[[1]][2, ])

  # no discounting in first period
  expect_equal(out$dcost[[1]][1, 1], 2)

  # discounting
  expect_equal(out$dcost[[1]][1, 2], 1/1.035)

  ## state health over time

  # same between simulations
  expect_equal(out$eff[[1]][1, ],
               out$eff[[1]][2, ])
  expect_equal(out$deff[[1]][1, ],
               out$deff[[1]][2, ])

  # no discounting in first period
  expect_equal(out$deff[[1]][1, 1], 1)

  # discounting
  expect_equal(out$deff[[1]][1, 2], 1/1.035)
})


test_that("discounting", {

  # no discounting
  out <-
    ce_sim(pop = pop,
           probs = probs,
           c_unit = c_unit(S, t_lim),
           e_unit = e_unit(S, t_lim),
           c_init = c_init(),
           delta = 0)

  expect_equal(out$cost[[1]],
               out$dcost[[1]])
  expect_equal(out$eff[[1]],
               out$deff[[1]])

  # 100% discounting i.e. halving
  out <-
    ce_sim(pop = pop,
           probs = probs,
           c_unit = c_unit(S, t_lim),
           e_unit = e_unit(S, t_lim),
           c_init = c_init(),
           delta = 1)

  expect_equal(out$dcost[[1]][1, ],
               c(2, 2^-(1:9)))
  expect_equal(out$deff[[1]][1, ],
               c(1, 2^-(1:9)))
})


test_that("unit values", {

})

#####################
# actual use case

test_that("only cost initial risk score", {

  # res_new <-
  #   ce_sim(probs,
  #          rc_unit(J,
  #                  c_icd_appt = 0,
  #                  c_icd_repl = 0,
  #                  c_shock = 0),
  #          e_unit,
  #          rc_init(c_icd = 0,
  #                  c_rscore = 1,
  #                  c_compl = 0),
  #          pdecr)
  #
  # c <-
  #   res_new$cost %>%
  #   map_dfc(rowSums) %>%
  #   as.matrix()
  #
  # # for 538 given icd
  # c_expected <-
  #   matrix(c(0, 538), ncol = 2, nrow = 1000, byrow = TRUE)
  #
  # expect_equal(c, c_expected, ignore_attr = TRUE)
  #
  # # check sign correct c1 - c0
  # m <- bcea(e, c, ref = 2)
  #
  # expect_equal(m$delta_c[[1]][1], 538 - 0)
})


##TODO:
# delta_c is difference between strategies in initial pop in ICD state
# init_risk6 - init_obs

# if all state occupancy utilities are 1 then total QALYs is t_lim*n_pop

