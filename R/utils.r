## General manipulation functions
make_non_neg <- function(x) {
  # makes a list of matrices non-negative
  #' x: list of matrices
  return(lapply(x, make_non_neg_inner))
}

make_non_neg_inner <- function(matrix) {
  # makes a matrix non-negative
  return(apply(matrix, 2, function(x) x + abs(min(0, min(x)))))
}

# function which takes a list L and vector v as input
# returns sum(v_i*L[[i]])
#' mat_list: list of k matrices
#' vec: a vector of length k
#' Output: a matrix of same dimensions as the entries in mat_list
star_prod <- function(vec, mat_list) {
  vec_mat <- 0
  for (i in seq_len(length(vec))) {
    if (vec[i] != 0) {
      vec_mat <- vec_mat + vec[i] * mat_list[[i]]
    }
  }
  return(vec_mat)
}

matrix_normalisation <- function(matrix) {
  # normalises a matrix so that the l1 norm of each column is 1
  normaliser <- diag(colSums(matrix))
  # solve Q
  normalised_matrix <- matrix %*% solve(normaliser)
  return(list(
    "normaliser" = normaliser,
    "normalised_matrix" = normalised_matrix
  ))
}

### functions for calculating JSD
jsd_calc <- function(x1, x2) {
  # calculate the Jensen Shannon divergence between
  # vectors x1 and x2
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
