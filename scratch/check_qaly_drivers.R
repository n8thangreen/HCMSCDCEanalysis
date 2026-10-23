suppressMessages(devtools::load_all("."))
load("data/jags_output_all_strat_direct.RData")
l <- mm1$BUGSoutput$sims.list$lambda
print(dim(l))
nm <- c("icd","risk_over_4","risk_over_6","no_icd")
try(load("data/transdat.RData")); print(names(transdat))
for (j in seq_along(names(transdat))) { cat(names(transdat)[j], "\n"); print(round(apply(l[, j, , ], c(2,3), mean), 4)) }
inp <- data.frame(t_repl=10,q_hcm=.88,u_icd=.9,u_shock=.875,u_implant=-.048,u_compl=-.096,p_compl=.043)
print(rowMeans(re_unit(30, inp, "icd"))[c(1,2,11,22,23,24)])
