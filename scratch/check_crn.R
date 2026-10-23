suppressMessages(devtools::load_all("."))
rpoint <- function(n, x) x
d <- readr::read_csv("inst/extdata/input_distns.csv", show_col_types = FALSE)
tv <- d[d$label=="tvicd",]; sc <- d[d$label=="sicd",]
set.seed(42); a <- do.call(rbind, lapply(1:1000, function(i) sample_param(tv)))
set.seed(42); b <- do.call(rbind, lapply(1:1000, function(i) sample_param(sc)))
cat("mean q_hcm  tv/sicd:", mean(a$q_hcm), mean(b$q_hcm), "\n")
cat("mean u_icd  tv/sicd:", mean(a$u_icd), mean(b$u_icd), "\n")
cat("mean u_shock tv/sicd:", mean(a$u_shock), mean(b$u_shock), "\n")
cat("cor(q_hcm) between devices:", cor(a$q_hcm, b$q_hcm), " first-iteration equal:", a$q_hcm[1]==b$q_hcm[1], "\n")
