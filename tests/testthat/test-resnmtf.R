data <- list(
  matrix(rnorm(100), 10, 10),
  matrix(rnorm(100), 10, 10)
)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


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
