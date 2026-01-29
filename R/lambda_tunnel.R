
#' convert BUGS transition probability matrix to
#' to model with tunnel state
#'
#' @param lambda S by S matrix
#' @param alpha proportion move from hcm to icd
#'
lambda_tunnel <- function(lambda, n_tunnel, alpha = 0) {
  # icd and shock state are duplicated
  n_states <- n_tunnel*2 + 4
  trans_mat <- matrix(0, ncol = n_states, nrow = n_states)

  # icd
  diag(trans_mat[, 2:n_tunnel]) <- lambda[1,1]
  trans_mat[n_tunnel, 1] <- lambda[1,1]
  pdeath <- lambda[1,3]

  # shock
  diag(trans_mat[, (n_tunnel + 1):(2*n_tunnel)]) <- lambda[1,2]
  # from shock to icd
  diag(trans_mat[(n_tunnel + 1):(2*n_tunnel), 3:n_tunnel]) <- 1 - pdeath
  trans_mat[2*n_tunnel - 1, 1] <- 1 - pdeath
  trans_mat[2*n_tunnel, 2] <- 1 - pdeath

  # all-cause mortality
  trans_mat[1:(2*n_tunnel), (2*n_tunnel) + 1] <- pdeath

  # sink states
  trans_mat[(2*n_tunnel) + 1, 2*n_tunnel + 1] <- 1
  trans_mat[(2*n_tunnel) + 3, 2*n_tunnel + 3] <- 1
  trans_mat[(2*n_tunnel) + 4, 2*n_tunnel + 4] <- 1

  # hcm
  trans_mat[(2*n_tunnel) + 2, 1] <- alpha*lambda[4,4]
  trans_mat[(2*n_tunnel) + 2, 2*n_tunnel + 2] <- (1 - alpha)*lambda[4,4]
  trans_mat[(2*n_tunnel) + 2, 2*n_tunnel + 3] <- lambda[4,5]
  trans_mat[(2*n_tunnel) + 2, 2*n_tunnel + 4] <- lambda[4,6]

  trans_mat
}

