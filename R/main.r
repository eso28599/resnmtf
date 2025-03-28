# ---------------------------------
# functions to run ResNMTF
# ---------------------------------


#' Calculate error
calculate_error <- function(data, current_f, current_s, current_g, n_v) {
  err <- c()
  for (v in 1:n_v) {
    x_hat <- current_f[[v]] %*% current_s[[v]] %*% t(current_g[[v]])
    err <- c(err, sum((data[[v]] - x_hat)**2) / sum((data[[v]])**2))
  }
  return(err)
}

normalisation_check <- function(current_f, current_g, current_s, n_v) {
  for (v in 1:n_v) {
    normal_f <- matrix_normalisation(current_f[[v]])
    current_f[[v]] <- normal_f$normalised_matrix
    normal_g <- matrix_normalisation(current_g[[v]])
    current_g[[v]] <- normal_g$normalised_matrix
    current_s[[v]] <- (normal_f$normaliser) %*%
      current_s[[v]] %*% normal_g$normaliser
  }
  return(list(
    "current_f" = current_f, "current_g" = current_g,
    "current_s" = current_s
  ))
}

# Application of ResNMTF to data!
#' @title Apply ResNMTF
#' @description Apply ResNMTF to data, without stability analysis and
#'              for a specific number of biclusters
#' @param data list of matrices, data to be factorised,
#' @param init_f list of matrices, initialisation for F matrices,
#' @param init_s list of matrices, initialisation for S matrices,
#' @param init_g list of matrices, initialisation for G matrices,
#' @param k_vec integer, vector of integers, number of clusters to consider,
#' @param phi list of matrices, default is NULL, restriction matrices for F,
#' @param xi list of matrices, default is NULL, restriction matrices for S,
#' @param psi list of matrices, default is NULL, restriction matrices for G,
#' @param n_iters integer, default is NULL, number of iterations to run for,
#' @param repeats integer, default is 5, minimum value of 2,
#'                number of repeats to use for ??
#' @param distance string, default is "euclidean",
#'                 distance metric to use within the bisilhouette score
#' @param no_clusts boolean, default is FALSE,
#'                  whether to return the factorisation
#'                  rather than the biclustering
#' @return list of results from ResNMTF
#' @export
res_nmtf_inner <- function(
    data,
    init_f = NULL,
    init_s = NULL,
    init_g = NULL, k_vec = NULL,
    phi = NULL, xi = NULL, psi = NULL,
    n_iters = NULL,
    repeats = 5,
    distance = "euclidean",
    no_clusts = FALSE) {
  n_v <- length(data)
  # initialise F, S and G based on svd decomposition if not given
  if (is.null(init_f) || is.null(init_g) || is.null(init_s)) {
    inits <- init_mats(data, k_vec)
    current_f <- inits$init_f
    current_s <- inits$init_s
    current_g <- inits$init_g
    current_lam <- inits$init_lambda
    current_mu <- inits$init_mu
  } else {
    # Take init_f, init_s, init_g as the initialised latent representations
    current_f <- init_f
    current_s <- init_s
    current_g <- init_g
    current_lam <- lapply(current_f, colSums)
    current_mu <- lapply(current_g, colSums)
  }
  # Update until convergence, or for n_iters times
  if (is.null(n_iters)) {
    total_err <- c()
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while ((err_diff > 1.0e-6)) {
      new_parameters <- update_matrices(
        x = data,
        input_f = current_f,
        input_s = current_s,
        input_g = current_g,
        lambda = current_lam,
        mu = current_mu,
        phi = phi,
        xi = xi,
        psi = psi
      )
      current_f <- new_parameters$output_f
      current_s <- new_parameters$output_s
      current_g <- new_parameters$output_g
      current_lam <- new_parameters$output_lam
      current_mu <- new_parameters$output_mu
      err <- calculate_error(data, current_f, current_s, current_g, n_v)
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err - err_temp)
      err_temp <- utils::tail(total_err, n = 1)
    }
  } else {
    total_err <- numeric(length = n_iters)
    for (t in 1:n_iters) {
      err <- numeric(length = length(current_f))
      new_parameters <- update_matrices(
        x = data,
        input_f = current_f,
        input_s = current_s,
        input_g = current_g,
        lambda = current_lam,
        mu = current_mu,
        phi = phi,
        xi = xi,
        psi = psi
      )
      current_f <- new_parameters$output_f
      current_s <- new_parameters$output_s
      current_g <- new_parameters$output_g
      current_lam <- new_parameters$output_lam
      current_mu <- new_parameters$output_mu
      total_err[t] <- mean(calculate_error(
        data, current_f,
        current_s, current_g, n_v
      ))
    }
  }
  normalised <- normalisation_check(current_f, current_g, current_s, n_v)
  current_f <- normalised[["current_f"]]
  current_g <- normalised[["current_g"]]
  current_s <- normalised[["current_s"]]
  # if only need to obtain factorisation, return values now
  if (no_clusts) {
    return(list(
      "output_f" = current_f, "output_s" = current_s,
      "output_g" = current_g
    ))
  }
  # find clustering results and bisilhouette score
  clusters <- obtain_biclusters(
    data, current_f,
    current_g, current_s, repeats, distance
  )
  if (is.null(n_iters)) {
    error <- mean(utils::tail(total_err, n = 10))
  } else {
    error <- utils::tail(total_err, n = 1)
  }
  return(list(
    "output_f" = current_f, "output_s" = current_s,
    "output_g" = current_g, "Error" = error,
    "All_Error" = total_err, "bisil" = clusters$bisil,
    "row_clusters" = clusters$row_clustering,
    "col_clusters" = clusters$col_clustering,
    "lambda" = current_lam,
    "mu" = current_mu
  ))
}

extract_bisils <- function(res_list, k_vec) {
  # extract scores
  err_list <- c()
  for (i in seq_along(k_vec)) {
    err_list <- c(err_list, res_list[[i]][["bisil"]])
  }
  return(err_list)
}

#' @param k_max integer, default is 6, must be greater than 2,
#'              largest value of k to be considered initially,
#' @param k_min integer, default is 3, must be greater than 1,
#'              smallest value of k to be considered initially,
#' @param distance string, default is "euclidean",
#'                 distance metric to use within the bisilhouette score
#' @param repeats integer, default is 5, minimum value of 2,
#'                number of repeats to use for ??
#' @param no_clusts boolean, default is FALSE, whether to return
#'                  only the factorisation or not,
#' @param sample_rate numeric, default is 0.9,
#'                    proportion of data to sample for stability analysis,
#' @param n_stability integer, default is 5,
#'                    number of times to repeat stability analysis,
#' @param stability boolean, default is TRUE,
#'                  whether to perform stability analysis or not,
#' @param stab_thres numeric, default is 0.4, threshold for stability analysis,
#' @param data list of matrices, data to be factorised,
#' @param init_f list of matrices, initialisation for F matrices,
#' @param init_s list of matrices, initialisation for S matrices,
#' @param init_g list of matrices, initialisation for G matrices,
#' @param k_vec integer, vector of integers, number of clusters to consider,
#' @param phi list of matrices, default is NULL, restriction matrices for F,
#' @param xi list of matrices, default is NULL, restriction matrices for S,
#' @param psi list of matrices, default is NULL, restriction matrices for G,
#' @param n_iters integer, default is NULL, number of iterations to run for,
#' @return list of results from ResNMTF
#' @export
#' @examples
#' data <- list(matrix(rnorm(100), nrow = 10), matrix(rnorm(100), nrow = 10))
#' apply_resnmtf(data = data, k_vec = c(3, 3), n_iters = 100)
#' apply_resnmtf(
#'   data = data, k_vec = c(3, 3), n_iters = 100,
#'   stability = FALSE
#' )
#' apply_resnmtf(data = data, k_vec = c(3, 3), n_iters = 100, no_clusts = TRUE)
#' apply_resnmtf(
#'   data = data, k_vec = c(3, 3),
#'   n_iters = 100, k_min = 3, k_max = 8
#' )
#' apply_resnmtf(
#'   data = data, k_vec = c(3, 3),
#'   n_iters = 100, k_min = 3, k_max = 8, distance = "euclidean"
#' )
#' apply_resnmtf(
#'   data = data, k_vec = c(3, 3),
#'   n_iters = 100, k_min = 3, k_max = 8,
#'   distance = "euclidean", repeats = 5
#' )
#' apply_resnmtf(
#'   data = data, k_vec = c(3, 3),
#'   n_iters = 100, k_min = 3, k_max = 8,
#'   distance = "euclidean", repeats = 5, no_clusts = FALSE
#' )
apply_resnmtf <- function(data, init_f = NULL, init_s = NULL,
                          init_g = NULL, k_vec = NULL,
                          phi = NULL, xi = NULL, psi = NULL,
                          n_iters = NULL, k_min = 3, k_max = 8,
                          distance = "euclidean", repeats = 5,
                          no_clusts = FALSE,
                          sample_rate = 0.9, n_stability = 5,
                          stability = TRUE, stab_thres = 0.4,
                          remove_unstable = TRUE) {
  # initialise phi etc matrices as zeros if not specified
  # otherwise multiply by given parameter
  n_v <- length(data)
  data <- make_non_neg(data)
  if (!typeof(data[[1]]) == "double") {
    data <- lapply(data, function(x) as.matrix(x))
  }
  # Normalise data
  data <- lapply(data, function(x) matrix_normalisation(x)$normalised_matrix)
  # initialise restriction matrices if not specified
  # views with no restrictions require no input
  phi <- init_rest_mats(phi, n_v)
  psi <- init_rest_mats(psi, n_v)
  xi <- init_rest_mats(xi, n_v)
  # if number of clusters has been specified method can be applied straight away
  if ((!is.null(k_vec))) {
    results <- res_nmtf_inner(
      data, init_f, init_s, init_g,
      k_vec, phi, xi, psi, n_iters,
      repeats, distance, no_clusts
    )
    # if using the original data, we want to perform stability analysis
    # otherwise we want the results
    if (stability) {
      return(stability_check(
        data, init_s, results,
        k_vec, phi, xi, psi, n_iters,
        repeats, no_clusts, distance, sample_rate,
        n_stability, stab_thres
      ))
    } else {
      return(results)
    }
  }
  # define set of k_s to consider
  k_vec <- k_min:k_max
  n_k <- length(k_vec)
  ones_vec <- rep(1, n_v)
  # apply method for each k to be considered
  # if on windows or linux operating system - do normal for loop
  if ((.Platform$OS.type == "windows") || (.Platform$OS.type == "unix")) {
    res_list <- vector("list", length = n_k)
    for (i in 1:n_k) {
      res_list[[i]] <- res_nmtf_inner(
        data, init_f, init_s, init_g,
        k_vec[i] * ones_vec, phi, xi, psi, n_iters,
        repeats, distance, no_clusts
      )
    }
  } else {
    # Register all the cores
    doParallel::registerDoParallel(min(parallel::detectCores(), length(k_vec)))
    res_list <- foreach::foreach(i = seq_along(k_vec)) %dopar% {
      res_nmtf_inner(
        data, init_f, init_s, init_g,
        k_vec[i] * ones_vec, phi, xi, psi, n_iters,
        repeats, distance, no_clusts
      )
    }
  }
  # extract scores
  err_list <- extract_bisils(res_list, k_vec)
  test <- k_vec[which.max(err_list)]
  max_k <- k_max
  # if best performing k is the largest k considered
  # apply method to k + 1 until this is no longer the case
  if (k_min != k_max) {
    while (test == max_k) {
      max_k <- max_k + 1
      k_vec <- c(k_vec, max_k)
      k <- max_k * rep(1, n_v)
      new_l <- length(k_vec)
      res_list[[new_l]] <- res_nmtf_inner(
        data, init_f, init_s, init_g,
        k, phi, xi, psi, n_iters,
        repeats, distance, no_clusts
      )
      err_list <- c(err_list, res_list[[new_l]][["bisil"]][1])
      test <- k_vec[which.max(err_list)]
    }
  }
  k <- which.max(err_list)
  results <- res_list[[k]]
  if (stability) {
    return(stability_check(
      data, init_s, results,
      k_vec[k], phi, xi, psi, n_iters,
      repeats, no_clusts, distance,
      sample_rate, n_stability,
      stab_thres, remove_unstable
    ))
  } else {
    return(results)
  }
}
