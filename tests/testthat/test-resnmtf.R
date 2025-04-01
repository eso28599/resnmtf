test_that("test_cond identifies rows/columns with only zeros", {
  # Case 1: Matrix with no zeros
  data <- list(matrix(1:9, 3, 3))
  expect_false(test_cond(data, attempt = 2))

  # Case 2: Matrix with a row of zeros
  data <- list(matrix(c(1, 2, 3, 0, 0, 0, 7, 8, 9), 3, 3, byrow = TRUE))
  expect_true(test_cond(data, attempt = 2))

  # Case 3: Matrix with a column of zeros
  data <- list(matrix(c(1, 0, 3, 4, 0, 6, 7, 0, 9), 3, 3, byrow = TRUE))
  expect_true(test_cond(data, attempt = 2))

  # Case 4: Matrix with all zeros
  data <- list(matrix(0, 3, 3))
  expect_true(test_cond(data, attempt = 2))

  # Case 5: Large matrix with no zeros
  data <- list(matrix(runif(100, 1, 10), 10, 10))
  expect_false(test_cond(data, attempt = 2))

  # Case 5: Large matrix with no zeros, but first attempt
  data <- list(matrix(runif(100, 1, 10), 10, 10))
  expect_true(test_cond(data, attempt = 1))
})

data <- list(
  matrix(rnorm(100), 10, 10),
  matrix(rnorm(100), 10, 10)
)

test_that("view shuffling works", {
  expect_equal(length(lapply(data, shuffle_view)), 2)
  expect_equal(dim(lapply(data, shuffle_view)[[1]])[1], 10)
})

data <- list(
  matrix(rnorm(2000), 100, 20),
  matrix(rnorm(2000), 100, 20)
)

row_clusters <- cbind(
  rbinom(100, 1, 0.5),
  rbinom(100, 1, 0.5),
  rbinom(100, 1, 0.5)
)
col_clusters <- cbind(
  rbinom(50, 1, 0.4),
  rbinom(50, 1, 0.4),
  rbinom(50, 1, 0.4)
)
n_col <- 50
data <- list(
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters) +
    abs(matrix(rnorm(100 * n_col), 100, n_col)),
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters) +
    abs(0.01 * matrix(rnorm(100 * n_col), 100, n_col))
)
test_that("get warning for negative matrix", {
  expect_warning(
    apply_resnmtf(-data[[1]], k_max = 4),
    "Matrix is not non-negative. Has been made non-negative."
  )
})

test_that("resnmtf runs", {
  results <- apply_resnmtf(data, k_max = 4)
  expect_equal(length(results$output_f), 2)
  expect_equal(relevance_results(
    results$row_clusters[[1]],
    results$col_clusters[[1]],
    row_clusters, col_clusters
  ), TRUE)
  relevance_results(
    results$row_clusters[[2]],
    results$col_clusters[[2]],
    row_clusters, col_clusters
  )
  expect_equal(relevance_results(
    results$row_clusters[[2]],
    results$col_clusters[[2]],
    row_clusters, col_clusters
  ), TRUE)
  expect_equal(dim(results$output_f[[1]])[1], 100)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_equal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_equal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_equal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_equal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
  print(sum(results$col_clusters[[2]][, 1] == col_clusters[, 1]))
  print(sum(results$col_clusters[[2]][, 1] == col_clusters[, 2]))
  print(sum(results$col_clusters[[2]][, 1] == col_clusters[, 3]))
})
