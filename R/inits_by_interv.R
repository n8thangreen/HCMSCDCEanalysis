
#' Start state populations
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(reshape2)
#'
#' data("ipd_risk")
#' inits_by_interv(
#'   ipd_risk,
#'   group_names = c("icd", "risk_over_6", "risk_over_4"))
#'
inits_by_interv <- function(ipd_risk,
                            group_names) {
  out <- ipd_risk

  # relabel logical to string
  for (i in group_names) {
    out <-
      mutate(
        out, "{i}" :=
          ifelse(!!sym(i),
                 yes = "ICD",
                 no = "low_risk"))
  }

  # group sizes
  out %>%
    select(group_names) %>%
    melt(id.vars = c(),
         variable.name = "risk",
         value.name = "status") %>%
    group_by(risk) %>%
    count(status)
}

