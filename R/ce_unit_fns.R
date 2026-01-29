
#' Sample costs
#' with tunnel states
#'
#' @param J max time (year)
#' @param inp unit values
#' @param name_interv
#'
#' @returns list
#' @export
rc_unit <- function(J, inp, name_interv) {

  c_icd_appt <- inp$c_appt
  c_icd_repl <- inp$c_icd_repl
  c_shock <- inp$c_shock
  t_repl <- inp$t_repl

  c_unit_score <- c(rep(2*c_icd_appt, t_repl),
                    rep(c_shock, t_repl), 0,
                    0, 0, 0)

  matrix(c_unit_score,
         ncol = J,
         nrow = length(c_unit_score),
         byrow = FALSE)
}


#' Sample starting costs
#' with tunnel states
#' @return list
#' @export
rc_init <- function(inp, name_interv) {

  c_compl <- c_complication(inp)
  c_icd <- inp$c_icd
  c_rscore <- inp$c_rscore
  c_compl <- inp$c_compl
  p_compl <- inp$p_compl
  c_icd_tot <- c_icd + p_compl*c_compl
  t_repl <- inp$t_repl
  n_tunnel <- t_repl - 1

  out <-
    switch(name_interv,
           icd =   c(c_icd_tot,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           risk_over_6 = c(c_icd_tot + c_rscore,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           risk_over_4 = c(c_icd_tot + c_rscore,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           any_rf = c(c_icd_tot + c_rscore,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           over_1_rf = c(c_icd_tot + c_rscore,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           over_2_rf = c(c_icd_tot + c_rscore,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           over_3_rf = c(c_icd_tot + c_rscore,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           no_icd = c(c_icd_tot + c_rscore,
                      rep(0, n_tunnel),
                      rep(0, t_repl), 0,
                      0, 0, 0))
  out
}


#' Sample starting health
#' with tunnel states
#' @return list
#' @export
re_init <- function(inp, name_interv) {

  u_implant <- inp$u_implant
  p_compl <- inp$p_compl
  u_compl <- inp$u_compl
  t_repl <- inp$t_repl
  n_tunnel <- t_repl - 1

  e_implant <- u_implant + p_compl*u_compl

  out <-
    switch(name_interv,
           icd =   c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           risk_over_6 = c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           risk_over_4 = c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           any_rf = c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           over_1_rf = c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           over_2_rf = c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           over_3_rf = c(e_implant,
                     rep(0, n_tunnel),
                     rep(0, t_repl), 0,
                     0, 0, 0),
           no_icd = c(e_implant,
                      rep(0, n_tunnel),
                      rep(0, t_repl), 0,
                      0, 0, 0))
  out
}

#' sample state health
#' with tunnel states
#' @param J max time (year)
#' @param inp
#' @param name_interv
#' @returns list
#' @export
re_unit <- function(J, inp, name_interv) {

  state_utils <-
    c(rep(inp$q_hcm*inp$u_icd, inp$t_repl),             # icd states
      rep(inp$q_hcm*inp$u_icd*inp$u_shock, inp$t_repl), # shock states
      0, inp$q_hcm, 0, 0)

  matrix(state_utils,
         ncol = J,
         nrow = length(state_utils),
         byrow = FALSE)
}


#' weighted infection, dislodgement
#' @export
c_complication <- function(inp) {

  ifelse(
    inp$p_inf_init + inp$p_dis_init == 0, 0,
    inp$p_inf_init*inp$c_inf + inp$p_dis_init*inp$c_dis)/(inp$p_inf_init + inp$p_dis_init)
}


#' sample from generic distribution
#' @export
sample_param <- function(inp) {
  sample_fns <- glue("r{inp$distn}(1, {inp$params})")
  rparam <- map(sample_fns, ~eval(parse(text = .)))
  names(rparam) <- inp$variable
  as.data.frame(rparam)
}

