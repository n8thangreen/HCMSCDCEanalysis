d <- read.csv("inst/extdata/scenario_input_values.csv", colClasses = "character", check.names = FALSE)
d$p_repl_compl <- ifelse(d$Scenario_Label == "tvicd", "0.0091", "0.002")
write.csv(d, "inst/extdata/scenario_input_values.csv", row.names = FALSE, quote = FALSE)
f <- read.csv("tests/testthat/scenario_input_values.csv", colClasses = "character", check.names = FALSE)
f$p_inf_yr <- "0"; f$p_dis_yr <- "0"; f$gamma_dis <- "1"
write.csv(f, "tests/testthat/scenario_input_values.csv", row.names = FALSE, quote = FALSE)
