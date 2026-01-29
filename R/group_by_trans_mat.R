
#' Aggregate transition matrix
#'
#' Create aggregate observed transition counts
#' see p.185 Table 5.2 chapter 5 in BCEA book
#'
#' Assume 3 state model:
#' 1) healthy 2) scd_count 3) non_scd_count
#'
#' @param data new state arrivals at time step
#'             result of running `annual_trans_counts()`
#'             e.g. data = data_rule$`FALSE`
#' @param max_year
#' @export
#'
group_by_trans_mat <- function(data,
                               max_year = NA) {

  trans_mat <- list()
  state_names <- c("healthy", "scd_count", "non_scd_count")
  max_year <- nrow(data) - 1

  # initialise counts
  trans_mat[[1]] <-
    rbind(c(0,0,0),
          c(0,0,0),
          c(0,0,0)) %>%
    `colnames<-`(state_names) %>%
    `rownames<-`(state_names)

  for (i in seq_len(max_year)) {

    trans_mat[[i + 1]] <- trans_mat[[i]]

    # new transitions
    trans_mat[[i + 1]]["healthy", state_names] <-
      trans_mat[[i + 1]]["healthy", state_names] +
      unlist(data[i, state_names])

    trans_mat[[i + 1]]["scd_count", "scd_count"] <-
      trans_mat[[i + 1]]["scd_count", "scd_count"] +
      trans_mat[[i]]["healthy", "scd_count"]

    trans_mat[[i + 1]]["non_scd_count", "non_scd_count"] <-
      trans_mat[[i + 1]]["non_scd_count", "non_scd_count"] +
      trans_mat[[i]]["healthy", "non_scd_count"]
  }

  names(trans_mat) <- c(0:max_year)
  trans_mat_total <- trans_mat[[as.character(max_year)]]

  cbind(trans_mat_total,
        n = rowSums(trans_mat_total))
}


#' Transition counts for multistate model
#' @return list of arrays
#'
trans_counts <- function(rs_name,
                         data,
                         CYCLE = 1) {
  transdat <-
    data %>%
    mutate(risk_status = ifelse(!!sym(rs_name),
                                yes = "ICD",
                                no = "low_risk"),
           risk_status = as.factor(risk_status), ) %>%
    split(.$risk_status) %>%     # convert to list
    map(.f = annual_trans_counts,
        cycle_length = CYCLE) %>%
    map(.f = group_by_trans_mat)

  # include zero matrix when no events
  if (rs_name == "no_icd") {
    transdat$ICD <- transdat$low_risk
    transdat$ICD[,] <- 0
  }

  transdat
}

