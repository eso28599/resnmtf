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

# number of rows/columns, divided by 3
n_col <- 60
n_row <- 60
row_clusters <- matrix(0, 3 * n_row, 3)
col_clusters <- matrix(0, 3 * n_col, 3)
for (i in 1:3) {
  row_clusters[((i - 1) * n_row + 1):(i * n_row), i] <- 1
  col_clusters[((i - 1) * n_col + 1):(i * n_col), i] <- 1
}

data <- list(
  row_clusters %*% diag(c(10, 10, 10)) %*% t(col_clusters) +
    0.1 * abs(matrix(rnorm(9 * n_row * n_col), 3 * n_row, 3 * n_col)),
  row_clusters %*% diag(c(10, 10, 10)) %*% t(col_clusters) +
    0.1 * abs(matrix(rnorm(9 * n_row * n_col), 3 * n_row, 3 * n_col))
)
test_that("get warning for negative matrix", {
  expect_warning(
    apply_resnmtf(list(-data[[1]]), k_max = 4),
    "Matrix is not non-negative. Has been made non-negative."
  )
})

# -----------------------------------
# tests with k specified
# -----------------------------------

test_that("resnmtf runs without stability with spurious removal", {
  results <- apply_resnmtf(data, k_vec = c(3, 3), stability = FALSE)
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], n_row * 3)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})

test_that("resnmtf runs with stability and spurious removal", {
  results <- apply_resnmtf(data, k_vec = c(3, 3))
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], n_row * 3)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})

test_that("resnmtf runs with stability and no spurious removal", {
  results <- apply_resnmtf(data,
    k_vec = c(3, 3),
    spurious = FALSE
  )
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], n_row * 3)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})

test_that("resnmtf runs with no stability and no spurious removal", {
  results <- apply_resnmtf(data,
    k_vec = c(3, 3),
    spurious = FALSE, stability = FALSE
  )
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], n_row * 3)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})
# -----------------------------------
# tests with k not specified
# -----------------------------------
test_that("resnmtf runs with k not specified", {
  results <- apply_resnmtf(data,
    k_max = 5,
    spurious = FALSE, stability = FALSE
  )
  expect_equal(length(results$output_f), 2)
  expect_equal(dim(results$output_f[[1]])[1], n_row * 3)
  expect_equal(dim(results$output_f[[1]])[2], 3)
  expect_setequal(colSums(results$row_clusters[[2]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[2]]), colSums(col_clusters))
  expect_setequal(colSums(results$row_clusters[[1]]), colSums(row_clusters))
  expect_setequal(colSums(results$col_clusters[[1]]), colSums(col_clusters))
})

# -----------------------------------
# tests for data reordering functions
# -----------------------------------
# generate output
data <- list(
  abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
  abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
)
row_names <- list(paste0("row_", 1:10), paste0("row_", 5:14))
col_names <- list(paste0("col_", 1:10), paste0("col_", 5:14))
n_views <- 2
for (view in 1:n_views) {
  rownames(data[[view]]) <- row_names[[view]]
  colnames(data[[view]]) <- col_names[[view]]
}
reordered_list <- reorder_data(data, n_views, row_names, col_names)

# perform tests
test_that("Row partition implemented correctly.", {
  expect_equal(length(Reduce(intersect, reordered_list$row_lists)), 0)
  expect_true(setequal(
    Reduce(union, reordered_list$row_lists),
    Reduce(union, row_names)
  ))
})
test_that("Column partition implemented correctly.", {
  expect_equal(length(Reduce(intersect, reordered_list$col_lists)), 0)
  expect_true(setequal(
    Reduce(union, reordered_list$col_lists),
    Reduce(union, col_names)
  ))
})
test_that("Row reordering implemented correctly.", {
  for (view in 1:n_views) {
    expect_setequal(
      rownames(data[[view]]),
      rownames(reordered_list$data_reordered[[view]])
    )
  }
})
test_that("Column reordering implemented correctly.", {
  for (view in 1:n_views) {
    expect_setequal(
      colnames(data[[view]]),
      colnames(reordered_list$data_reordered[[view]])
    )
  }
})

test_that("Returned to order correctly.", {
  n_col <- 60
  n_row <- 60
  row_clusters <- matrix(0, 3 * n_row, 3)
  col_clusters <- matrix(0, 3 * n_col, 3)
  for (i in 1:3) {
    row_clusters[((i - 1) * n_row + 1):(i * n_row), i] <- 1
    col_clusters[((i - 1) * n_col + 1):(i * n_col), i] <- 1
  }

  data <- list(
    row_clusters %*% diag(c(10, 10, 10)) %*% t(col_clusters) +
      0.1 * abs(matrix(rnorm(9 * n_row * n_col), 3 * n_row, 3 * n_col)),
    row_clusters %*% diag(c(10, 10, 10)) %*% t(col_clusters) +
      0.1 * abs(matrix(rnorm(9 * n_row * n_col), 3 * n_row, 3 * n_col))
  )
  row_names <- list(
    paste0("row_", 1:180),
    c(paste0("row_", 1:120), paste0("row_", 181:240))
  )
  col_names <- list(
    paste0("col_", 1:180),
    c(paste0("col_", 1:120), paste0("row_", 181:240))
  )
  for (view in 1:2) {
    rownames(data[[view]]) <- row_names[[view]]
    colnames(data[[view]]) <- col_names[[view]]
  }
  results <- apply_resnmtf(data,
    k_vec = c(3, 3),
    spurious = FALSE, stability = FALSE
  )
  for (i in 1:2) {
    expect_equal(rownames(results$output_f[[i]]), row_names[[i]])
    expect_equal(rownames(results$output_g[[i]]), col_names[[i]])
    expect_equal(rownames(results$row_clusters[[i]]), row_names[[i]])
    expect_equal(rownames(results$col_clusters[[i]]), col_names[[i]])
  }
})

# -----------------------------
# tests for naming functions
# -----------------------------
test_that(
  "Row names missing, row restriction provided but dimensions don't match.",
  {
    data <- list(
      abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
      abs(MASS::mvrnorm(20, mu = rep(0, 20), Sigma = diag(rep(1, 20))))
    )
    expect_error(
      give_names(data, 2, phi = matrix(1, 2, 2)),
      "Row restriction matrices implies shared rows between views
               with differing number of unnamed rows. Please name rows."
    )
  }
)
test_that("Row names missing but row restriction provided.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  row_names <- give_names(data, 2, phi = matrix(1, 2, 2))
  expect_equal(
    row_names$row_names[[1]], row_names$row_names[[2]]
  )
})

test_that(
  "Column names missing, column restriction provided, dimensions don't match",
  {
    data <- list(
      abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
      abs(MASS::mvrnorm(20, mu = rep(0, 20), Sigma = diag(rep(1, 20))))
    )
    expect_error(
      give_names(data, 2, psi = matrix(1, 2, 2)),
      "Column restriction matrices implies shared columns between
            views with differing number of unnamed columns.
            Please name columns."
    )
  }
)

test_that("Column names missing but column restriction provided.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  col_names <- give_names(data, 2, psi = matrix(1, 2, 2))
  expect_equal(
    col_names$col_names[[1]], col_names$col_names[[2]]
  )
})

test_that("One row missing a name correctly detected.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  rownames(data[[1]]) <- c(paste0("row_", 1:9), NA)
  rownames(data[[2]]) <- paste0("row_", 5:14)
  colnames(data[[1]]) <- paste0("col_", 1:10)
  colnames(data[[2]]) <- paste0("col_", 5:14)
  expect_error(
    give_names(data, 2),
    "Some rows missing names. Check row names."
  )
})
test_that("One view rows not named correctly detected.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  rownames(data[[1]]) <- paste0("row_", 1:10)
  colnames(data[[1]]) <- paste0("col_", 1:10)
  colnames(data[[2]]) <- paste0("col_", 5:14)
  expect_error(
    give_names(data, 2),
    "At least one view is missing row names. Please name missing rows."
  )
})
test_that("One column missing a name correctly detected.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  rownames(data[[1]]) <- paste0("row_", 1:10)
  rownames(data[[2]]) <- paste0("row_", 5:14)
  colnames(data[[1]]) <- c(paste0("col_", 1:9), NA)
  colnames(data[[2]]) <- paste0("col_", 5:14)
  expect_error(
    give_names(data, 2),
    "Some columns missing names. Check columns names."
  )
})
test_that("One view columns not named correctly detected.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  rownames(data[[1]]) <- paste0("row_", 1:10)
  rownames(data[[2]]) <- paste0("row_", 5:14)
  colnames(data[[1]]) <- paste0("col_", 1:10)
  expect_error(
    give_names(data, 2),
    "At least one view is missing column names. Please name missing columns."
  )
})
test_that("No row names present - names added.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  colnames(data[[1]]) <- paste0("col_", 1:10)
  colnames(data[[2]]) <- paste0("col_", 5:14)
  res <- give_names(data, 2)
  expect_equal(res$row_names[[1]], paste0("row_", 1:10))
  expect_equal(res$row_names[[2]], paste0("row_", 11:20))
})
test_that("No column names present - names added.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  rownames(data[[1]]) <- paste0("row_", 1:10)
  rownames(data[[2]]) <- paste0("row_", 5:14)
  res <- give_names(data, 2)
  expect_equal(res$col_names[[1]], paste0("col_", 1:10))
  expect_equal(res$col_names[[2]], paste0("col_", 11:20))
})

test_that("All rows/columns named - no change.", {
  data <- list(
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
    abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
  )
  rownames(data[[1]]) <- paste0("row_", 1:10)
  rownames(data[[2]]) <- paste0("row_", 5:14)
  colnames(data[[1]]) <- paste0("col_", 1:10)
  colnames(data[[2]]) <- paste0("col_", 5:14)
  res <- give_names(data, 2)
  expect_equal(res$row_names[[1]], paste0("row_", 1:10))
  expect_equal(res$row_names[[2]], paste0("row_", 5:14))
  expect_equal(res$col_names[[1]], paste0("col_", 1:10))
  expect_equal(res$col_names[[2]], paste0("col_", 5:14))
})
