# ---------------------------------
# utility functions
# ---------------------------------

#' @description makes a list of matrices non-negative
#' @param x list of matrices
#' @noRd
#' @return list of non-negative matrices
make_non_neg <- function(x) {
  return(lapply(x, make_non_neg_inner))
}

#' @description makes a matrix non-negative
#' @param matrix matrix to be made non-negative
#' @return non-negative matrix
#' @noRd
#' @details This function adds the absolute value of the minimum element
#'   of the matrix to each element of the matrix. This ensures that all
#'   elements of the matrix are non-negative.
make_non_neg_inner <- function(matrix) {
  # makes a matrix non-negative
  non_neg <- apply(matrix, 2, function(x) x + abs(min(0, min(x))))
  if (any(lapply(matrix, min) < 0)) {
    warning("Matrix is not non-negative. Has been made non-negative.")
  }
  return(non_neg)
}

#' @description Computes the product of a vector and a list of matrices
#' @param vec vector
#' @param mat_list list of matrices
#' @return matrix
#' @noRd
#' @details This function computes the product of a vector and a
#'          list of matrices where the product is defined as the sum
#'          of the product of each element of the vector
#'          with the corresponding matrix in the list. The function
#'          returns a matrix of the same dimensions as the entries in the list.
star_prod <- function(vec, mat_list) {
  vec_mat <- 0
  for (i in seq_len(length(vec))) {
    if (vec[i] != 0) {
      vec_mat <- vec_mat + vec[i] * mat_list[[i]]
    }
  }
  return(vec_mat)
}

#' @description normalises a matrix
#' @param matrix matrix to be normalised
#' @return list with two elements: normaliser and normalised_matrix
#' @noRd
#' @details This function normalises a matrix so that the l1 norm
#'          of each column is 1.
matrix_normalisation <- function(matrix) {
  normaliser <- diag(colSums(matrix))
  normalised_matrix <- matrix %*% solve(normaliser)
  return(list(
    "normaliser" = normaliser,
    "normalised_matrix" = normalised_matrix
  ))
}

#' @description Computes the Jensen-Shannon divergence between two vectors
#' @param x1 first vector
#' @param x2 second vector
#' @noRd
#' @return Jensen-Shannon divergence
jsd_calc <- function(x1, x2) {
  max_val <- max(x1, x2)
  d1 <- stats::density(x1, from = 0, to = max_val)
  d2 <- stats::density(x2, from = 0, to = max_val)
  d1$y[d1$x > max(x1)] <- 0
  d2$y[d2$x > max(x2)] <- 0
  return(
    suppressMessages(philentropy::JSD(rbind(d1$y, d2$y),
      unit = "log2", est.prob = "empirical"
    ))
  )
}

#' @description Computes the Jaccard similarity between two sets
#' @param a first set
#' @param b second set
#' @return Jaccard similarity
#' @noRd
#' @details This function computes the Jaccard similarity between
#'          two sets calculated as intersect(a,  b) / union(a, b).
#'          The function returns 0 if the union of the two sets is empty.
#'          The function returns 1 if the two sets are equal.
jaccard_func <- function(a, b) {
  # calculate jaccard between two vectors
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  if (union == 0) {
    return(0)
  } else {
    return(intersection / union)
  }
}

#' @description Computes the cartesian product of two sets
#' @param a first set
#' @param b second set
#' @noRd
#' @return a x b
cart_prod <- function(a, b) {
  # returns cartesian product of two sets
  prod <- c()
  # check a or b are not empty sets
  if (length(a) == 0 || length(b) == 0) {
    return(NULL)
  } else {
    for (k in seq_along(a)) {
      prod <- c(prod, paste(a[k], b))
    }
    return(prod)
  }
}

#' @title Calculate error
#' @description Calculate error for ResNMTF
#' @param data list of matrices, data to be factorised,
#' @param current_f list of matrices, current F matrices
#' @param current_s list of matrices, current S matrices
#' @param current_g list of matrices, current G matrices
#' @param n_v integer, number of views
#' @noRd
#' @return numeric vector, error for each view
calculate_error <- function(data, current_f, current_s, current_g, n_v) {
  err <- c()
  for (v in 1:n_v) {
    x_hat <- current_f[[v]] %*% current_s[[v]] %*% t(current_g[[v]])
    err <- c(err, sum((data[[v]] - x_hat)**2) / sum((data[[v]])**2))
  }
  return(err)
}

#' @title Normalise matrix
#' @description Normalises matrix factors
#' @param current_f list of matrices, current F matrices
#' @param current_g list of matrices, current G matrices
#' @param current_s list of matrices, current S matrices
#' @param n_v integer, number of views
#' @noRd
#' @return list of matrices, normalised F, G and S matrices
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


#' @title Extract bisilhouette scores
#' @description Extract bisilhouette scores from results list
#' @param res_list list of results from ResNMTF
#' @param k_vec vector of integers
#' @noRd
#' @return vector of bisilhouette scores
extract_bisils <- function(res_list, k_vec) {
  # extract scores
  err_list <- c()
  for (i in seq_along(k_vec)) {
    err_list <- c(err_list, res_list[[i]][["bisil"]])
  }
  return(err_list)
}

# -----------------------------------------------
# check functions
# -----------------------------------------------

#' @title Check if
#' @description Check if input is a positive integer
#' @param x numeric, input to check
#' @noRd
check_whole_number <- function(x, name) {
  if (!is.numeric(x)) {
    stop(paste(name, " must be a numeric."))
  }
  if (floor(x) != x || x <= 0) {
    stop(paste(name, " must be a positive integer."))
  }
}

#' @title Check integers
#' @description Check if inputs are positive integer
#' @param n_iters integer, number of iterations to run for
#' @param k_min integer, smallest value of k to be considered initially
#' @param k_max integer, largest value of k to be considered initially
#' @param num_repeats integer, number of num_repeats to use
#'                    for spurious bicluster removal
#' @param n_stability integer, number of times to repeat stability analysis
#' @noRd
check_integers <- function(n_iters, k_min, k_max,
                           num_repeats,
                           n_stability) {
  # check if integer inputs are integers
  if (!is.null(n_iters)) {
    check_whole_number(n_iters, "n_iters")
  }
  check_whole_number(num_repeats, "num_repeats")
  check_whole_number(n_stability, "n_stability")
  check_whole_number(k_min, "k_min")
  check_whole_number(k_max, "k_max")
  # check if k_max is greater than k_min
  if (k_max <= k_min) {
    stop("k_max must be greater than k_min.")
  }
}

#' @title Check boolean inputs
#' @description Check if inputs are boolean
#' @param no_clusts boolean, whether to return only the factorisation
#' @param stability boolean, whether to perform stability analysis or not
#' @param remove_unstable boolean, whether to remove unstable clusters or not
#' @noRd
check_boolean <- function(no_clusts,
                          stability, remove_unstable) {
  # check if boolean inputs are boolean
  if (!is.logical(no_clusts)) {
    stop("no_clusts must be a boolean.")
  }
  if (!is.logical(stability)) {
    stop("stability must be a boolean.")
  }
  if (!is.logical(remove_unstable)) {
    stop("remove_unstable must be a boolean.")
  }
}

#' @title Check numeric inputs
#' @description Check if inputs are numeric
#' @param sample_rate numeric, proportion of data to sample for
#'                    stability analysis
#' @param stab_thres numeric, threshold for stability analysis
#' @noRd
check_numeric <- function(sample_rate,
                          stab_thres) {
  # check if numeric inputs are numeric
  if (!is.numeric(sample_rate)) {
    stop("sample_rate must be a numeric.")
  }
  if (!is.numeric(stab_thres)) {
    stop("stab_thres must be a numeric.")
  }
  if (stab_thres < 0 || stab_thres > 1) {
    stop("stab_thres must be between 0 and 1.")
  }
  if (sample_rate <= 0 || sample_rate > 1) {
    stop("sample_rate must be greater than 0 and less than or equal to 1.")
  }
}

#' @title Check if inputs are list of matrices
#' @description Check if inputs are list of matrices
#' @param data list of matrices, data to be factorised
#' @param init_f list of matrices, initialisation for F matrices
#' @param init_s list of matrices, initialisation for S matrices
#' @param init_g list of matrices, initialisation for G matrices
#' @return list of matrices, data to be factorised
#' @noRd
check_lists <- function(data, init_f, init_s, init_g) {
  # check if data is a list of matrices
  if (!is.list(data)) {
    if (is.matrix(data)) {
      data <- list(data)
    } else {
      stop("Data must be a list of matrices or a matrix.")
    }
  }
  # check if init_f, init_s, init_g are lists of matrices
  if (!is.null(init_f)) {
    print(init_f)
    if (!is.list(init_f)) {
      stop("init_f must be a list of matrices or NULL.")
    }
  }
  if (!is.null(init_s)) {
    if (!is.list(init_s)) {
      stop("init_s must be a list of matrices or NULL.")
    }
  }
  if (!is.null(init_g)) {
    if (!is.list(init_g)) {
      stop("init_g must be a list of matrices or NULL.")
    }
  }
  return(data)
}

#' @title Check restriction matrix
#' @param matrix matrix, restriction matrix
#' @param name string, name of the matrix
#' @noRd
check_restriction_mat <- function(matrix, name) {
  if (!is.null(matrix)) {
    if (!is.matrix(matrix)) {
      stop(paste(name, "must be a matrix or NULL."))
    }
    if (any(matrix < 0)) {
      stop(paste(name, " must be a non-negative matrix."))
    }
    if (any(dim(matrix) != c(length(data), length(data)))) {
      stop(paste(name, " must be of the same dimensions as data."))
    }
  }
}

#' @title Check inputs
#' @description Check if inputs are valid,
#'              normalises and forces data to be non-negative
#' @param data list of matrices, data to be factorised
#' @param init_f list of matrices, initialisation for F matrices
#' @param init_s list of matrices, initialisation for S matrices
#' @param init_g list of matrices, initialisation for G matrices
#' @param k_vec vector of integers, number of clusters to consider
#' @param phi list of matrices, default is NULL, restriction matrices for F
#' @param xi list of matrices, default is NULL, restriction matrices for S
#' @param psi list of matrices, default is NULL, restriction matrices for G
#' @param n_iters integer, default is NULL, number of iterations to run for
#' @param k_min integer, smallest value of k to be considered initially
#' @param k_max integer, largest value of k to be considered initially
#' @param distance string, default is "euclidean",
#'                 distance metric to use within the bisilhouette score
#' @param num_repeats integer, default is 5, minimum value of 2,
#'                number of num_repeats to use for spurious bicluster removal
#' @param no_clusts boolean, default is FALSE, whether to return
#'                  only the factorisation or not
#' @param sample_rate numeric, default is 0.9,
#'                    proportion of data to sample for stability analysis
#' @param n_stability integer, default is 5,
#'                    number of times to repeat stability analysis
#' @param stability boolean, default is TRUE,
#'                  whether to perform stability analysis or not
#' @param stab_thres numeric, default is 0.4,
#'                   threshold for stability analysis
#' @param remove_unstable boolean, default is TRUE,
#'                        whether to remove unstable clusters or not
#' @return list of matrices, data to be factorised
#' @noRd
check_inputs <- function(data, init_f, init_s,
                         init_g, k_vec,
                         phi, xi, psi,
                         n_iters, k_min, k_max,
                         distance, num_repeats,
                         no_clusts,
                         sample_rate, n_stability,
                         stability, stab_thres,
                         remove_unstable) {
  # check lists
  data <- check_lists(data, init_f, init_s, init_g)
  # check integer inputs
  check_integers(
    n_iters, k_min, k_max,
    num_repeats, n_stability
  )
  # check boolean
  check_boolean(
    no_clusts, stability, remove_unstable
  )
  # check numeric
  check_numeric(
    sample_rate, stab_thres
  )
  # check if data is a non-negative matrix
  data <- make_non_neg(data)
  if (!typeof(data[[1]]) == "double") {
    warning("Data is not a double matrix. Converting to double.")
    data <- lapply(data, function(x) as.matrix(x))
  }
  # normalise matrix
  data <- lapply(data, function(x) matrix_normalisation(x)$normalised_matrix)
  ranks <- vapply(data, ncol, integer(1))
  # check if distance is a valid distance metric
  if (!distance %in% c("euclidean", "manhattan", "cosine")) {
    stop("distance must be one of 'euclidean', 'manhattan' or 'cosine'.")
  }
  # check if phi, xi, psi are matrices
  check_restriction_mat(phi, "phi")
  check_restriction_mat(xi, "xi")
  check_restriction_mat(psi, "psi")

  # check if k_vec is a vector of integers
  if (!is.null(k_vec)) {
    if (!is.numeric(k_vec)) {
      stop("k_vec must be a vector of integers.")
    }
    if (any(k_vec < 1)) {
      stop("k_vec must be a vector of integers greater than 1.")
    }
    if (length(k_vec) != length(data)) {
      stop("k_vec must be a vector of the same length as the number of views.")
    }
    if (any(k_vec > ranks)) {
      stop("k_vec must be a vector of integers less than or equal to the
            ranks of the views.")
    }
  } else {
    if (k_max > min(ranks)) {
      stop("k_max must be less than or equal to the minimum rank of the views.")
    }
  }
  return(data)
}
