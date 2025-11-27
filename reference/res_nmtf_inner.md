# ResNMTF inner functiion

Apply ResNMTF to data, without stability analysis and for a specific
number of biclusters

## Usage

``` r
res_nmtf_inner(
  data,
  row_indices,
  column_indices,
  init_f = NULL,
  init_s = NULL,
  init_g = NULL,
  k_vec = NULL,
  phi = NULL,
  xi = NULL,
  psi = NULL,
  n_iters = NULL,
  num_repeats = 5,
  spurious = TRUE,
  distance = "euclidean",
  no_clusts = FALSE
)
```

## Arguments

- data:

  list of matrices, data to be factorised

- row_indices:

  list of relevant row indices from all other views for each view

- column_indices:

  list of relevant column indices from all other views for each view

- init_f:

  list of matrices, initialisation for F matrices

- init_s:

  list of matrices, initialisation for S matrices

- init_g:

  list of matrices, initialisation for G matrices

- k_vec:

  integer, vector of integers, number of clusters to consider

- phi:

  list of matrices, default is NULL, restriction matrices for F

- xi:

  list of matrices, default is NULL, restriction matrices for S

- psi:

  list of matrices, default is NULL, restriction matrices for G

- n_iters:

  integer, default is NULL, number of iterations to run for

- num_repeats:

  integer, default is 5, minimum value of 2 number of num_repeats to use
  within spurious removal

- spurious:

  boolean, default is TRUE, whether or not spurious biclusters should be
  found and removed

- distance:

  string, default is "euclidean", distance metric to use within the
  bisilhouette score

- no_clusts:

  boolean, default is FALSE, whether to return the factorisation rather
  than the biclustering

## Value

list of results from ResNMTF
