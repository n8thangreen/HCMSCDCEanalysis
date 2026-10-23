suppressMessages(devtools::load_all("."))
v <- read.csv("inst/extdata/input_values.csv")
for (lbl in c("tvicd","sicd")) {
  g <- function(x) v$value[v$variable==x & v$label==lbl]
  inp <- list(p_inf_yr=g("p_inf_yr"), p_dis_yr=g("p_dis_yr"), gamma_dis=g("gamma_dis"), t_lead_cap=g("t_lead_cap"))
  pd <- p_lead_dis(inp, 30)
  cum <- function(n) 1 - prod((1-inp$p_inf_yr) * (1 - pd[1:n]))
  cat(lbl, "cumulative lead complications 7y:", round(cum(7),4), " 8y:", round(cum(8),4), " 10y:", round(cum(10),4), " 30y:", round(cum(30),4), "\n")
}
