# ---------------------------------
# functions to run ResNMTF
# ---------------------------------

#' @title ResNMTF inner functiion
#' @description Apply ResNMTF to data, without stability analysis and
#'              for a specific number of biclusters
#' @param data list of matrices, data to be factorised
#' @param row_indices list of relevant row indices from all other views
#'                    for each view
#' @param column_indices list of relevant column indices from all other views
#'    for each view
#' @param init_f list of matrices, initialisation for F matrices
#' @param init_s list of matrices, initialisation for S matrices
#' @param init_g list of matrices, initialisation for G matrices
#' @param k_vec integer, vector of integers, number of clusters to consider
#' @param phi list of matrices, default is NULL, restriction matrices for F
#' @param xi list of matrices, default is NULL, restriction matrices for S
#' @param psi list of matrices, default is NULL, restriction matrices for G
#' @param n_iters integer, default is NULL, number of iterations to run for
#' @param num_repeats integer, default is 5, minimum value of 2
#'                number of num_repeats to use within spurious removal
#' @param spurious boolean, default is TRUE, whether or not spurious biclusters
#'              should be found and removed
#' @param distance string, default is "euclidean",
#'                 distance metric to use within the bisilhouette score
#' @param no_clusts boolean, default is FALSE,
#'                  whether to return the factorisation
#'                  rather than the biclustering
#' @return list of results from ResNMTF
#' @export
res_nmtf_inner <- function(
    data, row_indices, column_indices,
    init_f = NULL, init_s = NULL, init_g = NULL,
    k_vec = NULL, phi = NULL, xi = NULL, psi = NULL,
    n_iters = NULL, num_repeats = 5, spurious = TRUE, distance = "euclidean",
    no_clusts = FALSE) {
  n_v <- length(data)
  initial_mats <- init_mats(
    data, n_v, k_vec,
    init_f, init_g, init_s
  )
  current_f <- initial_mats$current_f
  current_s <- initial_mats$current_s
  current_g <- initial_mats$current_g
  current_lam <- initial_mats$current_lam
  current_mu <- initial_mats$current_mu
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
        psi = psi,
        row_indices = row_indices,
        column_indices = column_indices
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
        psi = psi,
        row_indices = row_indices,
        column_indices = column_indices
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
    current_g, current_s, num_repeats, spurious, distance
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


#' @title Apply ResNMTF
#' @description Apply ResNMTF to data for a range of biclusters selecting
#'              the optimal number, with optional stability analysis
#' @param data list of n_v matrices, data to be factorised. If only one
#'             view is supplied, can be given as a matrix.
#' @param init_f list of matrices, initialisation for F matrices
#' @param init_s list of matrices, initialisation for S matrices
#' @param init_g list of matrices, initialisation for G matrices
#' @param k_vec vector of integers, number of clusters to consider
#'              in each view, default is NULL
#' @param phi n_v x n_v matrix, default is NULL, restriction matrices for F
#' @param xi n_v x n_v matrix, default is NULL, restriction matrices for S
#' @param psi n_v x n_v matrix, default is NULL, restriction matrices for G
#' @param n_iters integer, default is NULL, number of iterations to run for,
#'                otherwise will run until convergence
#' @param k_max positive integer, default is 6,
#'              largest value of k to be considered initially,
#' @param k_min positive integer, default is 3,
#'              smallest value of k to be considered initially,
#' @param distance string, default is "euclidean",
#'                 distance metric to use within the bisilhouette score
#' @param spurious boolean, default is TRUE, whether or not spurious biclusters
#'              should be found and removed
#' @param num_repeats integer, default is 5,
#'                number of repeats to use within stability analysis
#' @param no_clusts boolean, default is FALSE, whether to return
#'                  only the factorisation or not,
#' @param sample_rate numeric, default is 0.9,
#'                    proportion of data to sample for stability analysis,
#' @param n_stability integer, default is 5,
#'                    number of times to repeat stability analysis,
#' @param stability boolean, default is TRUE,
#'                  whether to perform stability analysis or not,
#' @param stab_thres numeric, default is 0.4, threshold for stability analysis,
#' @param remove_unstable boolean, default is TRUE,
#'                        whether to remove unstable clusters or not
#' @param use_parallel boolean, default is TRUE,
#'                     wheather to use parallelisation,
#'                     not applicable on Windows or linux machines
#' @return list of results from ResNMTF, containing the following:
#'         - output_f: list of matrices, F matrices
#'         - output_s: list of matrices, S matrices
#'         - output_g: list of matrices, G matrices
#'         - Error: numeric, mean error
#'         - All_Error: numeric, all errors
#'         - bisil: numeric, bisilhouette score
#'         - row_clusters: list of matrices, row clusters
#'         - col_clusters: list of matrices, column clusters
#'         - lambda: list of vectors, lambda vectors
#'         - mu: list of vectors, mu vectors
#' @export
#' @examples
#' row_clusters <- cbind(
#'   rbinom(100, 1, 0.5),
#'   rbinom(100, 1, 0.5),
#'   rbinom(100, 1, 0.5)
#' )
#' col_clusters <- cbind(
#'   rbinom(50, 1, 0.4),
#'   rbinom(50, 1, 0.4),
#'   rbinom(50, 1, 0.4)
#' )
#' n_col <- 50
#' data <- list(
#'   row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters) +
#'     abs(matrix(rnorm(100 * n_col), 100, n_col)),
#'   row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters) +
#'     abs(0.01 * matrix(rnorm(100 * n_col), 100, n_col))
#' )
#' apply_resnmtf(data, k_max = 4)
apply_resnmtf <- function(data, init_f = NULL, init_s = NULL,
                          init_g = NULL, k_vec = NULL,
                          phi = NULL, xi = NULL, psi = NULL,
                          n_iters = NULL, k_min = 3, k_max = 8,
                          distance = "euclidean", spurious = TRUE,
                          num_repeats = 5,
                          no_clusts = FALSE,
                          sample_rate = 0.9, n_stability = 5,
                          stability = TRUE, stab_thres = 0.4,
                          remove_unstable = TRUE, use_parallel = TRUE) {
  # initialise restriction matrices if not specified
  n_v <- length(data)
  # check naming
  named_data <- give_names(data, n_v, phi, psi)
  # get data indices
  reordering <- reorder_data(
    named_data$data, n_v, named_data$row_names, named_data$col_names
  )
  phi <- init_rest_mats(phi, n_v)
  psi <- init_rest_mats(psi, n_v)
  xi <- init_rest_mats(xi, n_v)
  # check inputs (normalises data)
  data <- check_inputs(
    named_data$data, init_f, init_s,
    init_g, k_vec,
    phi, xi, psi,
    n_iters, k_min, k_max,
    distance, num_repeats,
    no_clusts,
    sample_rate, n_stability,
    stability, stab_thres,
    remove_unstable, spurious
  )
  # if number of clusters has been specified method can be applied straight away
  if ((!is.null(k_vec))) {
    results <- res_nmtf_inner(
      data, reordering$row_indices, reordering$column_indices,
      init_f, init_s, init_g,
      k_vec, phi, xi, psi, n_iters,
      num_repeats, spurious, distance, no_clusts
    )
    # if using the original data, we want to perform stability analysis
    # otherwise we want the results
    if (stability) {
      results <- stability_check(
        data, results,
        k_vec, phi, xi, psi, n_iters,
        spurious, num_repeats,
        no_clusts, distance, sample_rate,
        n_stability, stab_thres
      )
    }
    return(results)
  }
  # define set of k_s to consider
  k_vec <- k_min:k_max
  n_k <- length(k_vec)
  ones_vec <- rep(1, n_v)
  # apply method for each k to be considered
  # if on windows or linux operating system - do normal for loop
  if (
    (.Platform$OS.type == "windows") || (.Platform$OS.type == "unix") ||
      (!use_parallel)
  ) {
    res_list <- vector("list", length = n_k)
    for (i in 1:n_k) {
      res_list[[i]] <- res_nmtf_inner(
        data, reordering$row_indices, reordering$column_indices,
        init_f, init_s, init_g,
        k_vec[i] * ones_vec, phi, xi, psi, n_iters,
        num_repeats, spurious, distance, no_clusts
      )
    }
  } else {
    # Register all the cores
    doParallel::registerDoParallel(min(parallel::detectCores(), length(k_vec)))
    res_list <- foreach::foreach(i = seq_along(k_vec)) %dopar% {
      res_nmtf_inner(
        data, reordering$row_indices, reordering$column_indices,
        init_f, init_s, init_g,
        k_vec[i] * ones_vec, phi, xi, psi, n_iters,
        num_repeats, spurious, distance, no_clusts
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
      new_len <- length(k_vec)
      res_list[[new_len]] <- res_nmtf_inner(
        data, reordering$row_indices, reordering$column_indices,
        init_f, init_s, init_g,
        k, phi, xi, psi, n_iters,
        num_repeats, spurious, distance, no_clusts
      )
      err_list <- c(err_list, res_list[[new_len]][["bisil"]][1])
      test <- k_vec[which.max(err_list)]
    }
  }
  k <- which.max(err_list)
  results <- res_list[[k]]
  if (stability) {
    results <- stability_check(
      data, results,
      k_vec[k], phi, xi, psi, n_iters,
      spurious, num_repeats,
      no_clusts, distance,
      sample_rate, n_stability,
      stab_thres, remove_unstable
    )
  }
  return(results)
}
