
##TODO: we need proper estimates of sd
##      email Menelous

# from: O’Mahony et al. Eur Heart J. 2014; 35(30)

hr <- list(
  age =
    c(mean = 0.98,
      sd = 1e-3), #(0.98 - 0.97)/2), ##TODO
  mwt =
    c(mean = 1.17,
      sd = 1e-3), #(1.17 - 1.01)/2),
  mwt2 =
    c(mean = 0.997,
      sd = 1e-3), #(0.997 - 0.99)/2),
  la =
    c(mean = 1.03,
      sd = 1e-3), #(1.03 - 1.01)/2),
  maxlvotg =
    c(mean = 1.004,
      sd = 1e-3), #(1.004 - 1.001)/2),
  fhxscd =
    c(mean = 1.58,
      sd = 1e-3), #(1.58 - 1.18)/2),
  nsvt =
    c(mean = 2.29,
      sd = 1e-3), #(2.29 - 1.64)/2),
  syncope =
    c(mean = 2.05,
      sd = 1e-3)) #(2.05 - 1.48)/2))


fixcoeff <- list(
  age =
    c(mean = -0.01799934,
      sd = 1e-3),
  mwt =
    c(mean = 0.15939858,
      sd = 1e-3),
  mwt2 =
    c(mean = -0.00294271,
      sd = 1e-3),
  la =
    c(mean = 0.0259082,
      sd = 1e-3),
  maxlvotg =
    c(mean = 0.00446131,
      sd = 1e-3),
  fhxscd =
    c(mean = 0.4583082,
      sd = 1e-3),
  nsvt =
    c(mean = 0.82639195,
      sd = 1e-3),
  syncope =
    c(mean = 0.71650361,
      sd = 1e-3))

