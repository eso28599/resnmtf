# ---------------------------------
# update functions
# ---------------------------------

#' @title initialise restriction matrices
#' @description if no restrictions required a matrix of zeros is returned,
#'              otherwise the restriction matrix is returned
#' @param mat matrix to be initialised
#' @param n_v number of views
#' @noRd
#' @return  restriction matrix, shape (n_v, n_v)
init_rest_mats <- function(mat, n_v) {
  # inititialise phi/psi/xi matrices
  # Check if the input matrix is NULL
  if (is.null(mat)) {
    # If NULL
    # initialize a new matrix of zeros with dimensions determined by n_v
    return(matrix(0, nrow = n_v, ncol = n_v))
  } else {
    # If not NULL, return the input matrix (validate the existing matrix)
    diag(mat) <- 0 # ensure diagonal is zero
    return(mat)
  }
}

#' @title initialise matrix factors
#' @description initialise matrix factors F, S, G etc.
#' @param x list of input data
#' @param n_v number of views
#' @param k_vec vector of cluster dimension
#' @param row_names list of row names for each view
#' @param col_names list of column names for each view
#' @noRd
#' @return a list of initial matrices, including
#'         inif_f, init_g, init_s, init_lambda, init_mu
init_mats <- function(
    x, n_v, k_vec,
    init_f, init_g, init_s) {
  row_names <- lapply(x, rownames)
  col_names <- lapply(x, colnames)
  # initialise F, S and G based on svd decomposition if not given
  if (is.null(init_f) || is.null(init_g) || is.null(init_s)) {
    inits <- init_mats_inner(x, k_vec, row_names, col_names, n_v)
    current_f <- inits$init_f
    current_s <- inits$init_s
    current_g <- inits$init_g
    current_lam <- inits$init_lambda
    current_mu <- inits$init_mu
  } else {
    # Take init_f, init_s, init_g as the initialised latent representations
    # check if they are valid
    current_f <- init_f
    current_s <- init_s
    current_g <- init_g
    current_lam <- lapply(current_f, colSums)
    current_mu <- lapply(current_g, colSums)
    for (i in 1:n_v) {
      rownames(current_f[[i]]) <- row_names[[i]]
      rownames(current_g[[i]]) <- col_names[[i]]
    }
  }
  return(list(
    "current_f" = current_f, "current_g" = current_g, "current_s" = current_s,
    "current_lam" = current_lam, "current_mu" = current_mu
  ))
}

#' @title inner function to initialise matrix factors
#' @description initialise matrix factors F, S, G etc.
#' @param x list of input data
#' @param k_vec vector of cluster dimension
#' @param sigma noise parameter, default is 0.05
#' @param row_names list of row names for each view
#' @param col_names list of column names for each view
#' @noRd
#' @return a list of initial matrices, including
#'         inif_f, init_g, init_s, init_lambda, init_mu
init_mats_inner <- function(
    x, k_vec, row_names, col_names, n_views,
    sigma = 0.05) {
  #' X: list of input data
  # Initialisation of F, S, G lists
  init_f <- vector("list", length = n_views)
  init_s <- vector("list", length = n_views)
  init_g <- vector("list", length = n_views)
  init_lambda <- vector("list", length = n_views)
  init_mu <- vector("list", length = n_views)
  # initialise based on svd
  for (i in 1:n_views) {
    k <- k_vec[i]
    vals <- 1:k
    ss <- svd(x[[i]])
    init_f[[i]] <- abs(ss$u[, vals])
    normal_f <- matrix_normalisation(init_f[[i]])
    init_f[[i]] <- normal_f$normalised_matrix
    init_s[[i]] <- abs(diag(ss$d)[vals, vals])
    init_s[[i]] <- init_s[[i]] + abs(MASS::mvrnorm(
      n = k,
      mu = rep(0, k), Sigma = sigma * diag(k)
    )[vals, vals])
    init_g[[i]] <- abs(ss$v[, vals])
    normal_g <- matrix_normalisation(init_g[[i]])
    init_g[[i]] <- normal_g$normalised_matrix
    init_s[[i]] <- (normal_f$normaliser) %*% init_s[[i]] %*% normal_g$normaliser
    init_lambda[[i]] <- colSums(init_f[[i]])
    init_mu[[i]] <- colSums(init_g[[i]])
    # add row and column names
    rownames(init_f[[i]]) <- row_names[[i]]
    rownames(init_g[[i]]) <- col_names[[i]]
  }

  return(list(
    "init_f" = init_f, "init_g" = init_g, "init_s" = init_s,
    "init_lambda" = init_lambda, "init_mu" = init_mu
  ))
}


#' @title Update F matrix
#' @description update F matrix based on the input data
#'              and the current matrix factors
#' @param x list of input data
#' @param input_f list of F matrices
#' @param input_s S^(v) matrix
#' @param input_g G^(v) matrix
#' @param lambda_in lambda^(v) vector
#' @param phi restriction matrix
#' @param v index of the view
#' @param row_indices list of relevant row indices from all other views
#' @noRd
#' @return updated F matrix
update_f <- function(
    x, input_f, input_s, input_g, lambda_in, phi,
    v, row_indices) {
  # Find numerator
  current_f <- input_f[[v]]
  numerator_matrix <- x %*% input_g %*% t(input_s)
  denominator_matrix <- current_f %*%
    input_s %*% t(input_g) %*% input_g %*% t(input_s)
  # update f
  phi_vec <- (phi + t(phi))[, v]
  lambda_mat <- 0.5 * t(matrix(lambda_in, length(lambda_in), nrow(current_f)))
  if (sum(phi_vec) == 0) {
    output_f <- current_f *
      ((numerator_matrix) / (denominator_matrix + lambda_mat))
  } else {
    num_mat_prod <- star_prod_relevant(phi_vec, input_f, current_f, row_indices)
    denom_mat_prod <- sum(phi_vec) * current_f
    output_f <- current_f * (
      (numerator_matrix + num_mat_prod) /
        (denominator_matrix + denom_mat_prod + lambda_mat)
    )
  }
  return(abs(output_f))
}

#' @title Update G matrix
#' @description update G matrix based on the input data
#'              and the current matrix factors
#' @param x list of input data
#' @param input_f F^(v) matrix
#' @param input_s S^(v) matrix
#' @param input_g list of G matrices
#' @param mu_in mu^(v) vector
#' @param psi restriction matrix
#' @param v index of the view
#' @param column_indices list of relevant column indices from all other views
#' @noRd
#' @return updated G matrix
update_g <- function(
    x, input_f, input_s, input_g, mu_in, psi,
    v, column_indices) {
  # Find numerator
  current_g <- input_g[[v]]
  numerator_matrix <- t(x) %*% input_f %*% input_s
  denominator_matrix <- current_g %*%
    t(input_s) %*% t(input_f) %*% input_f %*% input_s
  mu_mat <- 0.5 * t(matrix(mu_in, length(mu_in), nrow(current_g)))
  # update g
  if (sum(psi) == 0) {
    output_g <- current_g * (numerator_matrix / (denominator_matrix + mu_mat))
  } else {
    psi_vec <- (psi + t(psi))[, v]
    num_mat_prod <- star_prod_relevant(
      psi_vec, input_g, current_g,
      column_indices
    )
    denom_mat_prod <- sum(psi_vec) * current_g
    output_g <- current_g * (
      (numerator_matrix + num_mat_prod) /
        (denominator_matrix + denom_mat_prod + mu_mat)
    )
  }
  return(abs(output_g))
}

#' @title Update S matrix
#' @description update S matrix based on the input data
#'              and the current matrix factors
#' @param x list of input data
#' @param input_f F^(v) matrix
#' @param input_s list of S matrices
#' @param input_g G^(v) matrix
#' @param xi restriction matrix
#' @param v index of the view
#' @noRd
#' @return updated S matrix
update_s <- function(x, input_f, input_s, input_g, xi, v) {
  # Find numerator and demoninator
  current_s <- input_s[[v]]
  numerator_matrix <- t(input_f) %*% x %*% input_g
  denominator_matrix <- t(input_f) %*%
    input_f %*% current_s %*% t(input_g) %*% input_g
  # update s
  if (sum(xi) == 0) {
    output_s <- current_s * (numerator_matrix / denominator_matrix)
  } else {
    xi_vec <- (xi + t(xi))[, v]
    num_mat_prod <- star_prod(xi_vec, input_s)
    denom_mat_prod <- sum(xi_vec) * current_s
    output_s <- current_s * (
      (numerator_matrix + num_mat_prod) /
        (denominator_matrix + denom_mat_prod)
    )
  }
  return(abs(output_s))
}

#' @title Update lambda/mu vectors
#' @description update lambda/mu vectors based on the input data
#'              and the current matrix factors
#' @param vec vector to be updated
#' @param matrix matrix factor to be used for update
#' @noRd
#' @return updated vector
update_lm <- function(vec, matrix) {
  return(colSums(matrix) * vec)
}


#' @title Update matrices
#' @description update F, S, G matrices based on the input data
#'             and the current matrix factors
#' @param x list of input data
#' @param input_f list of F matrices
#' @param input_s list of S matrices
#' @param input_g list of G matrices
#' @param lambda lambda vector
#' @param mu mu vector
#' @param phi restriction matrix for F
#' @param xi restriction matrix for S
#' @param psi restriction matrix for G
#' @param row_indices list of relevant row indices from all other views
#'                    for each view
#' @param column_indices list of relevant column indices from all other views
#'    for each view
#' @noRd
#' @return list containing, ouput_f, output_s, output_g, output_lam, output_mu
update_matrices <- function(
    x, input_f, input_s, input_g, lambda, mu, phi, xi, psi,
    row_indices, column_indices) {
  n_v <- length(x)
  current_f <- input_f
  current_s <- input_s
  current_g <- input_g
  currentlam <- lambda
  currentmu <- mu
  # Update view-by-view
  for (v in 1:n_v) {
    # Update F
    current_f[[v]] <- update_f(
      x = x[[v]],
      input_f = current_f,
      input_s = current_s[[v]],
      input_g = current_g[[v]],
      lambda_in = currentlam[[v]],
      phi = phi,
      v = v,
      row_indices = row_indices[[v]]
    )
    # Update G
    current_g[[v]] <- update_g(
      x = x[[v]],
      input_f = current_f[[v]],
      input_s = current_s[[v]],
      input_g = current_g,
      mu_in = currentmu[[v]],
      psi = psi, v = v,
      column_indices = column_indices[[v]]
    )
    # Update S
    current_s[[v]] <- update_s(
      x = x[[v]],
      input_f = current_f[[v]],
      input_s = current_s,
      input_g = current_g[[v]],
      xi = xi, v = v
    )
    currentlam[[v]] <- update_lm(currentlam[[v]], current_f[[v]])
    currentmu[[v]] <- update_lm(currentmu[[v]], current_g[[v]])
  }
  return(list(
    "output_f" = current_f, "output_s" = current_s, "output_g" = current_g,
    "output_lam" = currentlam, "output_mu" = currentmu
  ))
}
