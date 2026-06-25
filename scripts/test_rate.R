load(here::here("data/jags_output.RData"))
print("Names in sims.list:")
print(names(mm1$BUGSoutput$sims.list))
if ("rate" %in% names(mm1$BUGSoutput$sims.list)) {
  print("Dim of rate:")
  print(dim(mm1$BUGSoutput$sims.list$rate))
  print("Summary of rate:")
  print(summary(mm1$BUGSoutput$sims.list$rate))
} else {
  print("rate is NOT in sims.list")
}
