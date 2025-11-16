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
  rbinom(200, 1, 0.5),
  rbinom(200, 1, 0.5),
  rbinom(200, 1, 0.5)
)
col_clusters <- cbind(
  rbinom(100, 1, 0.4),
  rbinom(100, 1, 0.4),
  rbinom(100, 1, 0.4)
)
n_col <- 100
data <- list(
  row_clusters %*% diag(c(10, 10, 10)) %*% t(col_clusters) +
    abs(matrix(rnorm(200 * n_col), 200, n_col)),
  row_clusters %*% diag(c(10, 10, 10)) %*% t(col_clusters) +
    abs(0.01 * matrix(rnorm(200 * n_col), 200, n_col))
)
test_that("get warning for negative matrix", {
  expect_warning(
    apply_resnmtf(list(-data[[1]]), k_max = 4),
    "Matrix is not non-negative. Has been made non-negative."
  )
})

test_that("resnmtf runs", {
  results <- apply_resnmtf(data, k_max = 3, stability = FALSE)
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], 200)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})

test_that("allowing no spurious removal runs ok", {
  results <- apply_resnmtf(data, k_max = 4, spurious = FALSE, stability = FALSE)
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], 200)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})
