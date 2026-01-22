# linkbrl

`linkbrl` provides functions for estimating **bipartite record linkage models** for categorical data using the **Classification Expectation Maximization (CEM)** algorithm. The package supports both the classical Fellegi–Sunter model and the graphical record linkage model of Steorts et al. (2016).

## Overview

The main function is `brl_cem`, which runs the CEM algorithm multiple times and outputs the best run. It can estimate the Fellegi-Sunter model (`model = "fs"`, with dedicated function `fs_brl_cem`)
and the graphical record linkage (`model = "graph"`, with dedicated function `graph_brl_cem`).

For both models, the *fast beta linkage (fabl)* prior on the matching structure is employed. At each iteration the CEM algorithm solves a linear optimization problem enforcing the at most one-to-one linking constraint of bipartite record linkage.

## Installation

The `linkbrl` R package is currently available only on GitHub. To install it, use the following code:

```r
# install.packages("remotes")
remotes::install_github("redivoe/linkbrl")
```

## Basic usage

Input datasets must be provided as integer-encoded categorical matrices, with the same variables (columns) in the same order.

```r
library(linkbrl)

fit <- brl_cem(X1, X2, model = "graph", reps = 5)

# Estimated coreference matrix
fit$Delta
# Predicted coreferent pairs
fit$coreferent_pairs
```

## References

- Kundinger, B., Reiter, J. P., and Steorts, R. C. (2025). *Efficient and Scalable Bipartite Matching with Fast Beta Linkage (fabl).* Bayesian Analysis, 20(3), 949–972. https://doi.org/10.1214/24-BA1427
- Steorts, R. C., Hall, R., and Fienberg, S. E. (2016). *A Bayesian Approach to Graphical Record Linkage and Deduplication.* Journal of the American Statistical Association, 111(516), 1660–1672. https://doi.org/10.1080/01621459.2015.1105807
- Redivo, E. (2026+). *Linking the Comparison and Graphical Approaches to Bipartite Matching.* Under review.
