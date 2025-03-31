# ---------------------------------
# utility functions
# ---------------------------------

#' @description makes a list of matrices non-negative
#' @param x list of matrices
#' @return list of non-negative matrices
make_non_neg <- function(x) {
  return(lapply(x, make_non_neg_inner))
}

#' @description makes a matrix non-negative
#' @param matrix matrix to be made non-negative
#' @return non-negative matrix
#' @details This function adds the absolute value of the minimum element
#'   of the matrix to each element of the matrix. This ensures that all
#'   elements of the matrix are non-negative.
make_non_neg_inner <- function(matrix) {
  # makes a matrix non-negative
  return(apply(matrix, 2, function(x) x + abs(min(0, min(x)))))
}

#' @description Computes the product of a vector and a list of matrices
#' @param vec vector
#' @param mat_list list of matrices
#' @return matrix
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
