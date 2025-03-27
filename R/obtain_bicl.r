# ---------------------------------
# functions to obtain biclusters
# and remove spurious ones
# ---------------------------------

#' @description shuffle a single view,
#'              ensuring returned matrix has no rows or column with all zeros
#' @param x_i matrix to shuffle
#' @return shuffled matrix
shuffle_view <- function(x_i) {
  dims <- dim(x_i)
  x_messed <- matrix(sample(x_i), dims[1], dims[2])
  while (
    any(colSums(x_messed) == 0) || any(rowSums(x_messed) == 0)
  ) {
    x_messed <- matrix(sample(x_i), dims[1], dims[2])
  }
  return(x_messed)
}

obtain_shuffled_f <- function(data, n_views, repeats, n_clusts) {
  f_mess <- vector(mode = "list", length = repeats)
  for (n in 1:repeats) {
    data_messed <- lapply(data, shuffle_view)
    f_mess[[n]] <- apply_resnmtf(
      data = data_messed,
      k_vec = n_clusts * rep(1, length = n_views),
      no_clusts = TRUE, stability = FALSE
    )$output_f
  }
  return(f_mess)
}

calculate_f_shuffle_jsd <- function(f_mess, i, j, repeats, n_clusts) {
  # calculate jsd between shuffled F matrices
  scores <- c()
  for (k in 1:n_clusts) {
    # jth repeat, ith view, kth cluster
    x1 <- ((f_mess[[j]])[[i]])[, k]
    for (l in (j + 1):repeats) {
      for (m in 1:n_clusts) {
        x2 <- ((f_mess[[l]])[[i]])[, m]
        scores <- c(scores, jsd_calc(x1, x2))
      }
    }
  }
  return(scores)
}

#' @title Get thresholds
#' @description obtain average and max JSD scores between
#'              shuffled F matrices for shuffled data
#' @param x list of x matrices
#' @param output_f list of F matrices
#' @param repeats number of repeats
#' @param n_views number of views
#' @param n_clusts number of biclusters
#' @return list containing average and max scores, as well as the shuffled x
get_thresholds <- function(x, output_f, repeats, n_views, n_clusts) {
  f_mess <- obtain_shuffled_f(x, n_views, repeats, n_clusts)
  avg_score <- c()
  max_score <- c()
  shuffled_f <- vector(mode = "list", length = n_views)
  for (i in 1:n_views) {
    scores <- c()
    for (j in 1:max(repeats - 1, 1)) {
      shuffled_f[[i]] <- cbind(shuffled_f[[i]], f_mess[[j]][[i]])
      scores <- c(scores, calculate_f_shuffle_jsd(
        f_mess, i, j, repeats, n_clusts
      ))
    }
    shuffled_f[[i]] <- cbind(shuffled_f[[i]], f_mess[[repeats]][[i]])
    avg_score <- c(avg_score, mean(scores))
    dens <- stats::density(scores)
    max_score <- c(max_score, dens$x[which.max(dens$y)])
  }
  return(list(
    "avg_score" = avg_score, "max_score" = max_score,
    "shuffled_f" = shuffled_f
  ))
}

check_biclusters <- function(data, output_f, repeats) {
  n_views <- length(data)
  n_clusts <- dim(output_f[[1]])[2]
  # updated results
  scores <- matrix(0, nrow = n_views, ncol = n_clusts)
  thresholds <- get_thresholds(data, output_f, repeats, n_views, n_clusts)
  for (i in 1:n_views) {
    x_noise <- thresholds$shuffled_f[[i]]
    for (k in 1:n_clusts) {
      x <- output_f[[i]][, k]
      scores[i, k] <- mean(apply(
        x_noise, 2,
        function(y) jsd_calc(x, y)
      ))
    }
  }
  return(list(
    "score" = scores,
    "avg_threshold" = thresholds$avg_score,
    "max_threshold" = thresholds$max_score
  ))
}

# Obtains biclustering from ResNMTF factorisation
# Args:
#  data: list of data matrices
#  output_f: list of F matrices
#  output_g: list of G matrices
#  output_s: list of S matrices
#  repeats: minimum value of 2
#  distance: distance metric to use, can be . Default is euclidean.
# Returns:
#  list containing row and column clusterings,
# and the bisilhouette score of the biclustering

obtain_biclusters <- function(data, output_f,
                              output_g, output_s, repeats, distance) {
  n_views <- length(output_f)
  row_clustering <- vector("list", length = n_views)
  col_clustering <- vector("list", length = n_views)

  # assign biclusters
  biclusts <- check_biclusters(data, output_f, repeats)
  for (i in 1:n_views) {
    row_clustering[[i]] <- apply(
      output_f[[i]],
      2, function(x) as.numeric(x > (1 / dim(output_f[[i]])[1]))
    )

    col_clustering[[i]] <- apply(
      output_g[[i]],
      2, function(x) as.numeric(x > (1 / dim(output_g[[i]])[1]))
    )
  }
  bisil <- c()
  # update realtions and
  # set biclusters that aren't strong enough to 0
  # and if bicluster is empty set row and cols to 0
  for (i in 1:n_views) {
    indices <- (
      ((biclusts$score[i, ]) < biclusts$max_threshold[i]) |
        ((biclusts$score[i, ]) == 0)
    )
    relations <- apply(output_s[[i]], 2, which.max)
    new_indices <- indices[relations] # i==0, col cluster i isn't a bicluster
    row_clustering[[i]] <- row_clustering[[i]][, relations]
    row_clustering[[i]][, new_indices] <- 0
    col_clustering[[i]][, new_indices] <- 0
    bisil <- c(
      bisil,
      bisilhouette::bisilhouette(
        data[[i]], row_clustering[[i]], col_clustering[[i]],
        method = distance
      )$bisil
    )
  }
  # calculate overall bisil
  bisil <- ifelse(sum(bisil) == 0, 0, mean(bisil[bisil != 0]))
  return(list(
    "row_clustering" = row_clustering,
    "col_clustering" = col_clustering, "bisil" = bisil
  ))
}
