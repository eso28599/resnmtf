# ResNMTF
Welcome to the landing page for the `resnmtf` package implementing Restrictive Non-Negative Matrix Tri-Factorisation from the paper "Multi-view biclustering via non-negative matrix tri-factorisation". 

# Installation 
You can install this R package using the following line of code:
```{r}
devtools::install_github("eso28599/bisilhouette") # requires devtools package to be installed
```

# Usage 

## Toy data
We define toy data with 2 views and 3 biclusters to apply ResNMTF to. We assume the biclusters have shared rows across views, but not shared columns.
```{r}
row_clusters <- cbind(
  rbinom(100, 1, 0.5),
  rbinom(100, 1, 0.5),
  rbinom(100, 1, 0.5)
)
col_clusters_v1 <- cbind(
  rbinom(50, 1, 0.4),
  rbinom(50, 1, 0.4),
  rbinom(50, 1, 0.4)
)
col_clusters_v2 <- cbind(
  rbinom(30, 1, 0.4),
  rbinom(30, 1, 0.4),
  rbinom(30, 1, 0.4)
)
data <- list(
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters_v2) +
    abs(matrix(rnorm(100 * 50), 100, 50)),
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters_v2) +
    abs(0.01 * matrix(rnorm(100 * 30), 100, 30))
)
```

The `apply_resnmtf` function . The data $X =\{X^{(1)}, \dots,X^{(n_v)}\}$ should be inputted as a list. 

The desired restrictions are imposed via restriction matrices $\phi$, $\psi$ and $\xi$ which enforce sharing row clusters, column clusters and row-column matching respectively. 
```{r}
phi <- matrix(0, 2, 2)
phi[1, 2] <- 1 # regularises towards shared rows between view 1 and view 2 
phi_val <- 200 # hyperparameter chosen to enforce regularisation
results <- apply_resnmtf(data, phi = phi_val * phi)
```
The multi-view biclustering is contained within two lists defining the row and column clusters associated with the biclusters. For a specific view, both are represented via binary logic matrices, with the $(i,j)^{th}$ element equal to $1$ if row/column $i$ belongs to bicluster $j$, and $0$ otherwise. 
```{r}
results$row_clusters
results$col_clusters
```

Note: non-zero $\xi$ restriction matrix also enforces the assumption that the scale of the biclusters is equal across views. 

### Single-view data
The method can also be applied on single-view data, inputted either as a list containing the matrix of the single-view:
```{r}
data <- list(
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters_v2) +
    abs(matrix(rnorm(100 * 50), 100, 50))
    )
results <- apply_resnmtf(data)
```
Or as a matrix:
```{r}
data <- row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters_v2) +
          abs(matrix(rnorm(100 * 50), 100, 50))
results <- apply_resnmtf(data)
```

# Citation
If you use our model in your work, please cite us with:

> Orme, E.S.C., Rodosthenous, T. and Evangelou, M., 2025. Multi-view biclustering via non-negative matrix tri-factorisation. arXiv preprint arXiv:2502.13698.

Or with the following bibtex entry:
```
@misc{resnmtf,
      title={Multi-view biclustering via non-negative matrix tri-factorisation}, 
      author={Ella S. C. Orme and Theodoulos Rodosthenous and Marina Evangelou},
      year={2025},
      eprint={2502.13698},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2502.13698}, 
}
```

