
#' Start state populations
#'
#' @param ipd_risk Dataframe containing patient risk data
#' @param group_names Vector of column names to group by
#' @importFrom dplyr mutate group_by count select summarise lag %>%
#' @importFrom rlang sym
#' @export
#'
#' @examples
#' library(dplyr)
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
  res_list <- lapply(group_names, function(i) {
    data.frame(
      risk = i,
      status = out[[i]]
    )
  })
  res_df <- do.call(rbind, res_list)

  res_df %>%
    group_by(risk) %>%
    count(status)
}

