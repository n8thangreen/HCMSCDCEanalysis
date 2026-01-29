
# convert from lognormal to normal parameters
# https://www.reddit.com/r/statistics/comments/8xy5qk/what_is_the_relationship_between_normal/
sd_norm <- function(x) unname(sqrt(log(x['sd']^2/(x['mean']^2) + 1)))
mean_norm <- function(x) unname(log(x['mean']) - (sd_norm(x)^2)/2)

# risk score
p_5_yr_scd <- function(dat, beta) {
  rf_names <- c("age","mwt","mwt2","la","maxlvotg",
                "fhxscd","nsvt","syncope")
  comp <- as.matrix(dat[, rf_names])%*%t(beta[, rf_names])
  1 - 0.998^exp(comp)
}

# from inv log scale
tranformLN <- function(hr)
  purrr::map(hr, ~c(mean = mean_norm(.), sd = sd_norm(.)))

# sample from frequentist distribution
rcoeff <- function(tranformLN, n_sample)
  purrr::map_df(tranformLN, ~rnorm(n_sample, .['mean'], .['sd']))


#
pop_sample_by_age <- function(risk_thresh,
                              data,
                              tmax = 2,       # time horizon
                              coef_samples) {
  p_risk <- list()

  # 5 year risk score for each sampled parameter set
  # and corresponding icd status
  # for increasing ages
  for (i in seq_len(n_sim)) {
    p_risk[[i]] <- list()
    rf <- data

    for (n in 1:tmax) {
      p_risk[[i]][[n]] <-
        p_5_yr_scd(rf, coeff_samples[i, ])

      # # point values from paper
      # p_5_yr_scd(covariates, map_dbl(fixcoeff, "mean"))

      is_icd <- p_risk[[i]][[n]] > risk_thresh
      rf <- rf[!is_icd, ]  # remove transitioned individuals
      rf$age <- rf$age + 1
    }
  }

  p_risk
}

