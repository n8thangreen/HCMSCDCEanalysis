d <- read.csv("inst/extdata/scenario_input_values.csv", colClasses = "character", check.names = FALSE)
names(d)[names(d) == "p_repl"] <- "p_repl_compl"
tv <- d$Scenario_Label == "tvicd"
d$p_inf_yr <- ifelse(tv, "0.002", "0.0005")
d$p_dis_yr <- ifelse(tv, "0.0029", "0.0027")
d$gamma_dis <- ifelse(tv, "1.6", "1")
write.csv(d, "inst/extdata/scenario_input_values.csv", row.names = FALSE, quote = FALSE)
