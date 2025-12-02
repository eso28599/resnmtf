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

#' @description Computes the product of a vector and a list of the
#'              relevant rows of matrices
#' @param vec vector
#' @param mat_list list of matrices
#' @param current_mat matrix of required dimension
#' @param indices list of indices to consider for each matrix in mat_list
#' @return matrix
#' @noRd
#' @details This function computes the product of a vector and a
#'          list of matrices where the product is defined as the sum
#'          of the product of each element of the vector
#'          with the relevant rows (as given by indices) of the corresponding
#'          matrix in the list. The function
#'          returns a matrix of the same dimensions as the entries in the list.
star_prod_relevant <- function(vec, mat_list, current_mat, indices) {
  vec_mat <- 0
  for (i in seq_len(length(vec))) {
    if (vec[i] != 0) {
      masked_matrix <- current_mat
      rownames(masked_matrix) <- rownames(current_mat)
      rows <- indices[[as.character(i)]]
      if (!any(is.na(rows))) {
        masked_matrix[rows, ] <- mat_list[[i]][rows, ]
        # vec_mat <- vec_mat + vec[i] * masked_matrix * length(rows)
        vec_mat <- vec_mat + vec[i] * masked_matrix * nrow(mat_list[[i]])
      }
    }
  }
  return(vec_mat / nrow(current_mat))
}

#' @description normalises a matrix
#' @param matrix matrix to be normalised
#' @return list with two elements: normaliser and normalised_matrix
#' @noRd
#' @details This function normalises a matrix so that the l1 norm
#'          of each column is 1.
matrix_normalisation <- function(matrix) {
  return(sweep(matrix, 2, colSums(matrix), FUN = "/"))
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
#' @param data_norms numeric vector, norms of each data matrix
#' @noRd
#' @return numeric vector, error for each view
calculate_error <- function(
    data, current_f, current_s, current_g, n_v, data_norms) {
  err <- rep(0, n_v)
  for (v in 1:n_v) {
    x_hat <- current_f[[v]] %*% current_s[[v]] %*% t(current_g[[v]])
    err[v] <- norm((data[[v]] - x_hat), "F")**2
  }
  err <- err / data_norms
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
    current_s[[v]] <- sweep(
      current_s[[v]], 2, colSums(current_f[[v]]) * colSums(current_g[[v]]),
      FUN = "*"
    )
    current_f[[v]] <- sweep(
      current_f[[v]], 2, colSums(current_f[[v]]),
      FUN = "/"
    )
    current_g[[v]] <- sweep(
      current_g[[v]], 2, colSums(current_g[[v]]),
      FUN = "/"
    )
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
#' @param spurious boolean, whether or not spurious biclusters should be
#'                 found and removed
#' @noRd
check_boolean <- function(no_clusts,
                          stability, remove_unstable, spurious) {
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
  if (!is.logical(spurious)) {
    stop("spurious must be a boolean.")
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
    if (is.list(init_f)) {
      stop("init_f must be a list of matrices or NULL.")
    }
  }
  if (!is.null(init_s)) {
    if (is.list(init_s)) {
      stop("init_s must be a list of matrices or NULL.")
    }
  }
  if (!is.null(init_g)) {
    if (is.list(init_g)) {
      stop("init_g must be a list of matrices or NULL.")
    }
  }
  return(data)
}

#' @title Check restriction matrix
#' @param matrix matrix, restriction matrix
#' @param name string, name of the matrix
#' @noRd
check_restriction_mat <- function(data, matrix, name) {
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
#' @param spurious boolean, whether or not spurious biclusters should be
#'                 found and removed
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
                         remove_unstable, spurious) {
  # check lists
  data <- check_lists(data, init_f, init_s, init_g)
  # check integer inputs
  check_integers(
    n_iters, k_min, k_max,
    num_repeats, n_stability
  )
  # check boolean
  check_boolean(
    no_clusts, stability, remove_unstable, spurious
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
  data <- lapply(data, function(x) matrix_normalisation(x))
  ranks <- vapply(data, ncol, integer(1))
  # check if distance is a valid distance metric
  if (!distance %in% c("euclidean", "manhattan", "cosine")) {
    stop("distance must be one of 'euclidean', 'manhattan' or 'cosine'.")
  }
  # check if phi, xi, psi are matrices
  check_restriction_mat(data, phi, "phi")
  check_restriction_mat(data, xi, "xi")
  check_restriction_mat(data, psi, "psi")

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

# -----------------------------------------------
# reordering and naming functions
# -----------------------------------------------

#' @title Give names to rows and columns
#' @description Give names to rows and columns if not provided,
#'              returns error if user input required
#' @param data list of matrices, data to be factorised
#' @param n_views integer, number of views
#' @param phi restriction matrix for F matrices, default is NULL
#' @param psi restriction matrix for G matrices, default is NULL
#' @return list with three elements: data, row_names, col_names
#' @noRd
give_names <- function(data, n_views, phi = NULL, psi = NULL) {
  names_missing_rows <- sapply(data, function(x) is.null(rownames(x)))
  names_missing_cols <- sapply(data, function(x) is.null(colnames(x)))
  # if all rows aren't named - name
  if (all(names_missing_rows)) {
    n <- 1
    for (i in 1:n_views) {
      if (is.null(rownames(data[[i]]))) {
        rownames(data[[i]]) <- paste0("row_", n:(n - 1 + nrow(data[[i]])))
        n <- n + nrow(data[[i]])
      }
      # if restriction matrices provided, check number of rows match and
      # assumes they are well aligned
      if (!is.null(phi)) {
        for (j in (min(i + 1, n_views)):n_views) {
          if (phi[i, j] > 0 && nrow(data[[i]]) != nrow(data[[j]])) {
            stop("Row restriction matrices implies shared rows between views
               with differing number of unnamed rows. Please name rows.")
          } else if (phi[i, j] > 0) {
            # if restriction matrix implies shared rows, copy names
            rownames(data[[j]]) <- rownames(data[[i]])
          }
        }
      }
    }
    # if some views are completely missing names, user must specify
  } else if (any(names_missing_rows)) {
    stop("At least one view is missing row names. Please name missing rows.")
  } else {
    # if some rows are missing names within a view but others are provided
    # user must specify
    if (any(sapply(data, function(x) any(is.na(rownames(x)))))) {
      stop("Some rows missing names. Check row names.")
    }
  }
  # now check columns
  if (all(names_missing_cols)) {
    n <- 1
    for (i in 1:n_views) {
      if (is.null(colnames(data[[i]]))) {
        colnames(data[[i]]) <- paste0("col_", n:(n - 1 + ncol(data[[i]])))
        n <- n + ncol(data[[i]])
      }
      # if restriction matrices provided, check number of rows match and
      # assumes they are well aligned
      if (!is.null(psi)) {
        for (j in (min(i + 1, n_views)):n_views) {
          if (psi[i, j] > 0 && ncol(data[[i]]) != ncol(data[[j]])) {
            stop("Column restriction matrices implies shared columns between
            views with differing number of unnamed columns.
            Please name columns.")
          } else if (psi[i, j] > 0) {
            # if restriction matrix implies shared columns, copy names
            colnames(data[[j]]) <- colnames(data[[i]])
          }
        }
      }
    }
  } else if (any(names_missing_cols)) {
    stop(
      "At least one view is missing column names. Please name missing columns."
    )
  } else {
    # if any columns are missing names
    if (any(sapply(data, function(x) any(is.na(colnames(x)))))) {
      stop("Some columns missing names. Check columns names.")
    }
  }
  return(list(
    "data" = data,
    "row_names" = lapply(data, rownames),
    "col_names" = lapply(data, colnames)
  ))
}

#' @title Produce indices of shared rows and columns
#' @description Produce indices of shared rows and columns
#' @param row_lists list of vectors, lists of row names for each existing
#'                  view subset
#' @param col_lists list of vectors, lists of column names for each existing
#'                  view subset
#' @param existing_view_subsets_r list of vectors, lists of existing
#'                                view subsets for rows
#' @param existing_view_subsets_c list of vectors, lists of existing
#'                                view subsets for columns
#' @param n_views integer, number of views
#' @return list with two elements: shared_rows, shared_cols.
#'         Each element is a list corresponding to a view where each entrie
#'         is a hash where the keys are the other views and
#'         the values are the shared row/column names between the two views.
#' @noRd
produce_indices <- function(
    row_lists, col_lists,
    existing_view_subsets_r, existing_view_subsets_c, n_views) {
  shared_rows <- vector("list", length = n_views)
  shared_cols <- vector("list", length = n_views)
  for (view1 in 1:n_views) {
    shared_v1 <- hash::hash()
    shared_c_v1 <- hash::hash()
    # add in n_views ==1 case
    for (view2 in (1:n_views)[-view1]) {
      common_rows <- unlist(
        row_lists[
          sapply(
            existing_view_subsets_r,
            function(x) all(c(view1, view2) %in% x)
          )
        ]
      )
      common_cols <- unlist(
        col_lists[
          sapply(
            existing_view_subsets_c,
            function(x) all(c(view1, view2) %in% x)
          )
        ]
      )
      if (length(common_rows) != 0) {
        shared_v1[[as.character(view2)]] <- common_rows
      } else {
        shared_v1[[as.character(view2)]] <- NA
      }
      if (length(common_cols) != 0) {
        shared_c_v1[[as.character(view2)]] <- common_cols
      } else {
        shared_c_v1[[as.character(view2)]] <- NA
      }
    }
    shared_rows[[view1]] <- shared_v1
    shared_cols[[view1]] <- shared_c_v1
  }
  return(list("shared_rows" = shared_rows, "shared_cols" = shared_cols))
}

#' @title Reorder data
#' @description Reorder data so that rows and columns belonging to the same
#'              view subsets are grouped together
#' @param data list of matrices, data to be factorised
#' @param n_views integer, number of views
#' @param row_names list of vectors, row names for each view
#' @param col_names list of vectors, column names for each view
#' @return list with two elements:
#'         row_indices: a list where each element is a hash corresponding to a
#'          view, where the keys are the other views and
#'          the values are the shared row names between the two views.
#'         col_indices:
#'           a list where each element is a hash corresponding to a
#'          view, where the keys are the other views and
#'          the values are the shared column names between the two views
#' @noRd
reorder_data <- function(data, n_views, row_names, col_names) {
  # for each view v, find R^(v) - the list of existing view subsets A for view v
  # define relevant power set, not including the empty subset
  # a row is in existing view subset A
  power_set <- rje::powerSetCond(1:n_views)
  existing_view_subsets_r <- list()
  existing_view_subsets_c <- list()
  row_lists <- list()
  col_lists <- list()
  # extract relevant row_names
  for (views in power_set) {
    # define {1, ..., n} \ A
    negate_v_sub <- (1:n_views)[!((1:n_views) %in% views)]
    # find rows belonging to views v in A
    rows <- Reduce(intersect, row_names[views])
    # find rows belonging to views v not in A
    rows_negation <- Reduce(union, row_names[negate_v_sub])
    # select the rows present in A but not in {1, ..., n} \ A
    rows_in_a <- setdiff(rows, rows_negation)
    if (length(rows_in_a) != 0) {
      existing_view_subsets_r <- c(existing_view_subsets_r, list(views))
      row_lists <- c(row_lists, list(rows_in_a))
    }
    # find cols belonging to views v in A
    cols <- Reduce(intersect, col_names[views])
    # find cols belonging to views v not in A
    cols_negation <- Reduce(union, col_names[negate_v_sub])
    # select the cols present in A but not in {1, ..., n} \ A
    cols_in_a <- setdiff(cols, cols_negation)
    if (length(cols_in_a) != 0) {
      existing_view_subsets_c <- c(existing_view_subsets_c, list(views))
      col_lists <- c(col_lists, list(cols_in_a))
    }
  }

  indices <- produce_indices(
    row_lists, col_lists,
    existing_view_subsets_r, existing_view_subsets_c, n_views
  )
  return(list(
    "row_indices" = indices$shared_rows,
    "col_indices" = indices$shared_cols
  ))
}
