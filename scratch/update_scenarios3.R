d <- read.csv("inst/extdata/scenario_input_values.csv", colClasses = "character", check.names = FALSE)
i <- match("u_compl", names(d))
d <- cbind(d[1:i], u_repl = "-0.024", d[(i+1):ncol(d)])
d$t_lead_cap <- "10"
write.csv(d, "inst/extdata/scenario_input_values.csv", row.names = FALSE, quote = FALSE)
