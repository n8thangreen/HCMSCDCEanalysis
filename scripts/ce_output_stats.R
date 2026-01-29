
# cost-effectiveness summary stats
# and tables

library(BCEA)

load(file = "data/bcea_data.RData")


# tables
BCEA:::ce_table(he)
BCEA:::summary.bcea(he)  ##TODO: error
BCEA::tabulate_means(he)


# proportions shocked

sum(res$pop$risk6[2, , 1])/filter(init, risk == "risk_over_6", status == "ICD")$n
# 0.229149
sum(res$pop$obs[2, , 1])/filter(init, risk == "icd", status == "ICD")$n
# 0.1866873
sum(res$pop$risk4[2, , 1])/filter(init, risk == "risk_over_4", status == "ICD")$n

