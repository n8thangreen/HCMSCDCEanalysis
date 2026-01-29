
# HCM SCD risk prediction model
# N Green, UCL
#
# prep individual level data to use in BUGS model
#
# the individual-level dataset is split into those who
# are given the intervention (ICD) and those who are not
# for different decision rules.
#
# the individual data are then summed into annual counts
# from the health state.
#
# then aggregated across some time horizon and formatted
# as a transition matrix


library(dplyr)
library(purrr)

## STATA data
# mydata <- haven::read_dta("raw data/hcmdata.dta")
# # single imputed data set
# data_set1 <- mydata %>% filter(set == 1)
# save(data_set1, file = "data/data_set1.RData")

data("data_set1")

CYCLE <- 1 #year

# non-deterministic decision rule
# random at boundary
FUZZY_RISK <- FALSE

fuzzy_noise <-
  if (FUZZY_RISK) {
    rnorm(nrow(data_set1), 0, 0.001)
  } else 0

## risk factors:
# mwt30  : max wall thickness
# nsvt   : NSVT
# fhxscd : family history of SCD
# syncope: unexplained syncope
# no noise on decision

# probability of 5 year risk >6%, >4%
ipd_risk <-
  data_set1 %>%
  mutate(risk_over_6 =
           data_set1$risk_5_years + fuzzy_noise > 0.06,
         risk_over_4 =
           data_set1$risk_5_years + fuzzy_noise > 0.04,
         num_rf = mwt30 + nsvt + fhxscd + syncope,
         no_icd = FALSE,
         any_rf = num_rf > 0,
         over_1_rf = num_rf > 1,
         over_2_rf = num_rf > 2,
         over_3_rf = num_rf > 3)

# complete data
save(ipd_risk, file = "data/ipd_risk.RData")


#############################################
# aggregate transition counts by risk group

is_icd_names <-
  c("risk_over_6", "risk_over_4", "any_rf", "over_1_rf",
    "over_2_rf", "over_3_rf", "icd", "no_icd")

transdat <- list()

for (i in is_icd_names) {
  transdat[[i]] <- trans_counts(i, ipd_risk)
}

library(reshape2)
trans_long <-  melt(transdat)

save(transdat, file = "data/transdat.RData")

