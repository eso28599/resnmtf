# ---------------------------------
# functions for stability selection
# ---------------------------------


#' @title calculate Jaccard index between two biclusterings
#' @description calculate Jaccard index between two biclusterins
#' @param row_c row clusters
#' @param col_c column clusters
#' @param true_r true row clusters
#' @param true_c true column clusters
#' @param m number of row clusters
#' @param n number of column clusters
#' @noRd
#' @return matrix of Jaccard indices
jaccard_main <- function(row_c, col_c, true_r, true_c, m, n) {
  # initialise storage of jaccard index between pairs
  samps <- seq_along(row_c[, 1])
  feats <- seq_along(col_c[, 1])
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    r_i <- samps[row_c[, i] == 1]
    c_i <- feats[col_c[, i] == 1]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n) {
      tr_i <- samps[true_r[, j] == 1]
      tc_i <- feats[true_c[, j] == 1]
      m_j <- cart_prod(tr_i, tc_i)
      jac_mat[i, j] <- jaccard_func(m_i, m_j)
    }
  }
  return(jac_mat)
}

#' @title calculate the relevance of biclusters
#' @description calculate the relevannce of biclusters obtained from
#'              sub-sampled data, using the initial biclustering
#'              as the ground truth
#' @param row_c row clusters
#' @param col_c column clusters
#' @param true_r true row clusters
#' @param true_c true column clusters
#' @noRd
#' @return vector of relevance values for each bicluster
relevance_results <- function(row_c, col_c, true_r, true_c) {
  m <- ncol(row_c)
  n <- ncol(true_r)
  m_0 <- sum(colSums(row_c) != 0) # no of clusters actually detected
  n_0 <- sum(colSums(true_r) != 0) # no of true clusters

  # edge cases
  # if either:
  #     - no biclusters detected but some are present
  #     - no biclusters present but some are detected
  # return 0
  if ((m_0 == 0 && n_0 != 0) || (n_0 == 0 && m_0 != 0)) {
    return(0)
  }
  # if no biclusters present and none detected
  # return 1
  if (m_0 == 0 && n_0 == 0) {
    return(1)
  }
  # initialise storage of jaccard index between pairs
  jac_mat <- jaccard_main(row_c, col_c, true_r, true_c, m, n)
  return(apply(jac_mat, 2, max))
}

#' Test whether there are any columns/rows with only zeros
#' @description Test whether there are any columns/rows with only zeros
#' @param data list of matrices to be tested
#' @param attempt number of attempts made to sample data
#' @noRd
#' @return TRUE if any columns/rows with only zeros, FALSE otherwise
test_cond <- function(data, attempt) {
  if (attempt == 1) {
    return(TRUE)
  }
  return(any(unlist(lapply(
    data,
    function(x) {
      any(colSums(x) == 0) | any(rowSums(x) == 0)
    }
  ))))
}

#' @title number of biclusters
#' @description determine number of biclusters
#' @param results results of apply_resnmtf/res_nmtf_inner
#' @noRd
#' @return number of biclusters
number_biclusters <- function(results) {
  return(sum(as.numeric(lapply(
    results$row_clusters,
    function(x) sum(colSums(x))
  ))))
}

#' @title shuffle view
#' @description shuffle view
#' @param data list of matrices to be shuffled
#' @param new_data list of shuffled matrices
#' @param dims dimensions of the view to be shuffled
#' @param dim dimensions of the first view
#' @param dim_1 dimensions of the first view
#' @param i index of the view to be shuffled
#' @param sample_rate rate at which to sample
#' @param row_samples list of rows to sample
#' @param col_samples list of cols to sample
#' @noRd
initial_shuffle <- function(
    data, new_data, dims, dim, i,
    sample_rate, row_samples, col_samples) {
  if ((dims[1]) == dim[1]) {
    row_samples[[i]] <- row_samples[[1]]
  } else {
    row_samples[[i]] <- sample(dims[1], (dims[1] * sample_rate))
  }
  if ((dims[2]) == dim[2]) {
    col_samples[[i]] <- col_samples[[1]]
  } else {
    col_samples[[i]] <- sample(dims[2], (dims[2] * sample_rate))
  }
  new_data[[i]] <- data[[i]][row_samples[[i]], col_samples[[i]]]
  rownames(new_data[[i]]) <- rownames(data[[i]])[row_samples[[i]]]
  colnames(new_data[[i]]) <- colnames(data[[i]])[col_samples[[i]]]
  return(list(
    "new_data" = new_data,
    "row_samples" = row_samples,
    "col_samples" = col_samples
  ))
}

#' @title check for empty matrix
#' @description check if matrix has any rows/columns with only zeros
#' @param data list of matrices to be checked
#' @param i index of the view to be checked
#' @noRd
#' @return TRUE if any rows/columns with only zeros, FALSE otherwise
check_empty <- function(data, i) {
  return(any(colSums(data[[i]]) == 0) || any(rowSums(data[[i]]) == 0))
}

#' @title sample view
#' @description sample view
#' @param data list of matrices to be shuffled
#' @param new_data list of shuffled matrices
#' @param dims dimensions of the view to be shuffled
#' @param dim_1 dimensions of the first view
#' @param i index of the view to be shuffled
#' @param row_samples list of rows to sample
#' @param col_samples list of cols to sample
#' @param sample_rate rate at which to sample
#' @noRd
#' @return list of shuffled matrices (new_data),
#'         rows (row_samples) and columns (col_samples)
sample_view <- function(data, i, new_data, dim_1,
                        row_samples, col_samples, sample_rate) {
  dims <- dim(data[[i]])
  shuffled <- initial_shuffle(
    data, new_data, dims, dim_1, i, sample_rate, row_samples, col_samples
  )
  new_data <- shuffled[["new_data"]]
  row_samples <- shuffled[["row_samples"]]
  col_samples <- shuffled[["col_samples"]]
  if (check_empty(new_data, i)) {
    zeros_cols <- colSums(new_data[[i]]) != 0
    zeros_rows <- rowSums(new_data[[i]]) != 0
    if ((dims[1]) == dim_1[1]) {
      for (p in 1:i) {
        row_samples[[p]] <- row_samples[[p]][zeros_rows]
      }
    } else {
      row_samples[[i]] <- row_samples[[i]][zeros_rows]
    }
    if ((dims[2]) == dim_1[2]) {
      for (p in 1:i) {
        col_samples[[p]] <- col_samples[[p]][zeros_cols]
      }
    } else {
      col_samples[[i]] <- col_samples[[i]][zeros_cols]
    }
    for (p in 1:i) {
      new_data[[p]] <- data[[p]][row_samples[[p]], col_samples[[p]]]
    }
  }
  return(list(
    "new_data" = new_data,
    "row_samples" = row_samples,
    "col_samples" = col_samples
  ))
}

#' @title stability repeat
#' @description perform stability analysis on the results
#' @param results results of apply_resnmtf/res_nmtf_inner
#' @param data list of matrices to be shuffled
#' @param dim_1 dimensions of the first view
#' @param k number of biclusters
#' @param phi phi restriction parameter
#' @param xi xi restriction parameter
#' @param psi psi restriction parameter
#' @param n_iters number of iterations
#' @param num_repeats number of repeats in removing spurious biclusters
#' @param distance distance metric to be used within res_nmtf_inner
#' @param spurious boolean, whether or not spurious biclusters
#'                 should be found and removed
#' @param n_views number of views
#' @param sample_rate rate at which to sample
#' @noRd
#' @return matrix (i, j) of relevance values
#'         for each bicluster (j) in each view (i)
#'         and whether stability analysis was performed
#'         (stability_performed)
stability_repeat <- function(results, data, dim_1, k, phi, xi, psi, n_iters,
                             num_repeats, distance, spurious,
                             n_views, sample_rate) {
  new_data <- vector(mode = "list", length = n_views)
  row_samples <- vector(mode = "list", length = n_views)
  col_samples <- vector(mode = "list", length = n_views)
  relevance <- matrix(0, nrow = n_views, ncol = k)
  # turn this into a function to be used with lapply
  attempt <- 1
  # need to check that
  while (test_cond(new_data, attempt)) {
    if (attempt == 20) {
      print("Unable to perform stability analysis due to sparsity of data.")
      return(list("stability_performed" = FALSE))
    }
    row_samples[[1]] <- sample(dim_1[1], (dim_1[1] * sample_rate))
    col_samples[[1]] <- sample(dim_1[2], (dim_1[2] * sample_rate))
    new_data[[1]] <- data[[1]][row_samples[[1]], col_samples[[1]]]
    if (check_empty(new_data, 1)) {
      zeros_cols <- colSums(new_data[[1]]) != 0
      zeros_rows <- rowSums(new_data[[1]]) != 0
      row_samples[[1]] <- row_samples[[1]][zeros_rows]
      col_samples[[1]] <- col_samples[[1]][zeros_cols]
      new_data[[1]] <- data[[1]][row_samples[[1]], col_samples[[1]]]
      rownames(new_data[[1]]) <- rownames(data[[1]])[row_samples[[1]]]
      colnames(new_data[[1]]) <- colnames(data[[1]])[col_samples[[1]]]
    }
    if (n_views > 1) {
      for (i in 2:n_views) {
        shuffled <- sample_view(
          data, i, new_data, dim_1, row_samples, col_samples, sample_rate
        )
        new_data <- shuffled[["new_data"]]
        row_samples <- shuffled[["row_samples"]]
        col_samples <- shuffled[["col_samples"]]
      }
    }
    attempt <- attempt + 1
  }
  new_results <- res_nmtf_inner(
    new_data,
    k_vec = k * rep(1, n_views),
    phi = phi, xi = xi, psi = psi,
    n_iters = n_iters, num_repeats = num_repeats,
    spurious = spurious,
    distance = distance
  )
  # extract results
  for (i in 1:n_views) {
    relevance[i, ] <- relevance[i, ] + relevance_results(
      new_results$row_clusters[[i]],
      new_results$col_clusters[[i]],
      results$row_clusters[[i]][row_samples[[i]], ],
      results$col_clusters[[i]][col_samples[[i]], ]
    )
  }
  return(list("relevance" = relevance, "stability_performed" = TRUE))
}

#' @title stability check
#' @description perform stability analysis on the biclustering results
#' @param data list of input matrices
#' @param results results of apply_resnmtf/res_nmtf_inner
#' @param k number of biclusters found
#' @param phi phi restriction parameter
#' @param xi xi restriction parameter
#' @param psi psi restriction parameter
#' @param n_iters number of iterations
#' @param spurious boolean, whether or not spurious biclusters
#'              should be found and removed
#' @param num_repeats number of repeats in removing spurious biclusters
#' @param no_clusts number of clusters
#' @param distance distance metric to be used within res_nmtf_inner
#' @param sample_rate rate at which to sample
#' @param n_stability number of stability checks to perform
#' @param stab_thres threshold for stability
#' @param remove_unstable whether to remove unstable biclusters
#' @noRd
#' @return results of biclustering with unstable biclusters removed
#'         (if applicable) or a list containing the initial results and
#'         relevance values for each bicluster
stability_check <- function(data, results,
                            k, phi, xi, psi, n_iters,
                            spurious, num_repeats,
                            no_clusts, distance, sample_rate = 0.9,
                            n_stability = 5, stab_thres = 0.6,
                            remove_unstable = TRUE) {
  # check whether stability check needs to be performed
  if (number_biclusters(results) == 0) {
    print("No biclusters detected!")
    return(results)
  }
  n_views <- length(data)
  dim_1 <- dim(data[[1]])
  # initialise storage of results
  relevance <- matrix(0, nrow = n_views, ncol = k)
  for (t in 1:n_stability) {
    stab_repeat <- stability_repeat(
      results, data, dim_1, k, phi, xi, psi, n_iters,
      num_repeats, distance, spurious, n_views, sample_rate
    )
    relevance <- relevance + stab_repeat$relevance
    if (!stab_repeat$stability_performed) {
      return(results)
    }
  }
  relevance <- relevance / n_stability
  if (!remove_unstable) {
    return(list("res" = results, "relevance" = relevance))
  } else {
    for (i in 1:n_views) {
      # set clusters not deemed stable to have 0 members
      results$row_clusters[[i]][, relevance[i, ] < stab_thres] <- 0
      results$col_clusters[[i]][, relevance[i, ] < stab_thres] <- 0
    }
    return(results)
  }
}
