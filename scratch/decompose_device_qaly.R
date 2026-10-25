suppressMessages(devtools::load_all("."))
load("data/jags_output_all_strat_direct.RData"); load("data/transdat.RData")
l <- mm1$BUGSoutput$sims.list$lambda; j <- which(names(transdat)=="icd")
lam <- apply(l[, j, , ], c(2,3), mean)
v <- read.csv("inst/extdata/input_values.csv")
run <- function(lbl, mod = list(), J = 30) {
  inp <- as.list(setNames(v$value[v$label==lbl], v$variable[v$label==lbl])); if (is.null(inp$u_compl)) inp$u_compl <- -0.096; inp[names(mod)] <- mod
  T <- inp$t_repl_max; S <- 2*T+4
  P <- array(lambda_tunnel(lam, T, 0), dim=c(S,S,1,1))
  pop <- array(0, dim=c(S,J,1)); pop[1,1,1] <- 1
  o <- ce_sim(list(icd=pop), list(icd=P), inp=as.data.frame(inp), c_unit=rc_unit, e_unit=re_unit, c_init=rc_init, e_init=re_init, distn=FALSE, interv_names="icd")
  c(QALY=o$e_total[1], cost=o$c_total[1])
}
vars <- list(base=list(), no_u_repl=list(u_repl=0), no_lead=list(p_inf_yr=0,p_dis_yr=0), no_repl_compl=list(p_repl_compl=0), same_t_repl=list(t_repl_max=10))
res <- sapply(vars, function(m) c(tvicd=run("tvicd", m)["QALY"], sicd=run("sicd", m)["QALY"]))
rownames(res) <- c("tvicd","sicd"); print(round(res,4)); print(round(res["sicd",]-res["tvicd",],4))
cat("time in state by device (alive-years):\n")
