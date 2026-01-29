
# raw data summary


library(dplyr)
library(purrr)


data("data_set1")

data_set1 <-
  data_set1 |>
  mutate(non_scd = as.numeric((cvs_all & !scd) | non_cvs),
         cens = as.numeric(!non_scd & !scd),
         composite = anyrip_tx | scd)


# follow up period
summary(data_set1$time)

# composite endpoint
sum(data_set1$scdtype != 0, na.rm = TRUE)
mean(data_set1$scdtype != 0, na.rm = TRUE)

# actual SCD
sum(data_set1$scdtype == 3, na.rm = TRUE)
mean(data_set1$scdtype == 3, na.rm = TRUE)

# appropriate shock
sum(data_set1$scdtype == 1, na.rm = TRUE)
mean(data_set1$scdtype == 1, na.rm = TRUE)

# survived cardiac arrest
sum(data_set1$scdtype == 2, na.rm = TRUE)
mean(data_set1$scdtype == 2, na.rm = TRUE)

# all cv deaths combined
sum(data_set1$cvs_all, na.rm = TRUE)
mean(data_set1$cvs_all, na.rm = TRUE)

sum(data_set1$cvd, na.rm = TRUE)
mean(data_set1$cvd, na.rm = TRUE)

sum(data_set1$cvs_other, na.rm = TRUE)
mean(data_set1$cvs_other, na.rm = TRUE)

# all deaths and transplants
sum(data_set1$anyrip_tx, na.rm = TRUE)
mean(data_set1$anyrip_tx, na.rm = TRUE)

# non-scd death
sum(data_set1$non_scd, na.rm = TRUE)
mean(data_set1$non_scd, na.rm = TRUE)

# non-cardiovascular death
sum(data_set1$non_cvs, na.rm = TRUE)
mean(data_set1$non_cvs, na.rm = TRUE)

sum(data_set1$cens, na.rm = TRUE)
mean(data_set1$cens, na.rm = TRUE)

sum(data_set1$composite, na.rm = TRUE)
mean(data_set1$composite, na.rm = TRUE)
