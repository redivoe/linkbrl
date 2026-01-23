#' Simulate data under a graphical record linkage model
#'
#' Simulate data for bipartite record linkage following the graphical record linkage model
#' of Steorts, Hall, and Fienberg (2016).
#'
#' @param n1 Integer. Number of units in dataset 1.
#' @param n2 Integer. Number of units in dataset 2.
#' @param n12 Integer. Number of linked units shared across datasets. These are
#'   assumed to correspond to rows \code{1:n12} in \code{X1} and \code{X2}.
#' @param p Integer. Number of variables (columns).
#' @param L Integer vector of length \code{p}. \code{L[h]} is the number of
#'   categories for variable \code{h}. Categories are encoded as integers
#'   \code{1, ..., L[h]}.
#' @param params Optional list with elements \code{theta} and \code{beta}. If
#'   provided, these parameters are used instead of being randomly generated.
#'   \code{theta} should have dimensions \code{c(p, max(L))} with only
#'   columns \code{1:L[h]} used for feature \code{h}. \code{beta} should be a
#'   numeric vector of length \code{p} with entries in \code{[0, 1]}.
#'
#' @return A list with components:
#' \describe{
#'   \item{X1}{Integer matrix \code{n1 x p} of categorical values.}
#'   \item{X2}{Integer matrix \code{n2 x p} of categorical values.}
#'   \item{beta}{Numeric vector of length \code{p} of distortion probabilities.}
#'   \item{theta}{Numeric array \code{c(p, max(L))} of category probabilities.}
#'   \item{L}{The input vector with the number of categories for each variable.}
#'   \item{Delta}{A sparse \code{TsparseMatrix} of dimension \code{n1 x n2} with
#'     ones on the first \code{n12} diagonal entries, encoding the assumed pairing
#'     between units \code{1:n12} across datasets.}
#' }
#'
#' @details
#' When \code{params} is missing, \code{theta[h, 1:L[h]]} is sampled from a
#' symmetric Dirichlet distribution and
#' \code{beta[h]} is sampled from a \code{Beta(2, 20)}.
#'
#' @references
#' Steorts, R. C., Hall, R., and Fienberg, S. E. (2016).
#' \emph{A Bayesian Approach to Graphical Record Linkage and Deduplication}.
#' Journal of the American Statistical Association, 111(516), 1660--1672.
#' \doi{10.1080/01621459.2015.1105807}
#'
#' @examples
#'   set.seed(1)
#'   data <- sim_graph_brl(n1 = 20, n2 = 15, n12 = 10, p = 3, L = c(2, 3, 4))
#'   str(data)
#'
#' @export
sim_graph_brl <- function(n1, n2, n12, p, L, params = NULL){

  stopifnot("sim_graph_brl(): n2 must be <= n1." = n2 <= n1)

  if(is.null(params)){
    theta <- array(dim = c(p, max(L)))
    for(h in 1:p){
      theta[h, 1:L[h]] <- c(MCMCpack::rdirichlet(n = 1, alpha = rep(1, L[h])))
    }
    beta <- rbeta(n = p, shape1 = 2, shape2 = 20)
  }else{
    theta <- params$theta
    beta <- params$beta
  }

  X1 <- array(dim = c(n1, p))
  X2 <- array(dim = c(n2, p))
  X12 <- array(dim = c(n12, p))
  for(h in 1:p){
    if(n12 < n1){
      X1[(n12+1):n1, h] <- sample.int(n = L[h], size = n1 - n12, replace = TRUE, prob = theta[h, 1:L[h]])
    }
    if(n12 < n2){
      X2[(n12+1):n2, h] <- sample.int(n = L[h], size = n2 - n12, replace = TRUE, prob = theta[h, 1:L[h]])
    }
    X12[, h] <- sample.int(n = L[h], size = n12, replace = TRUE, prob = theta[h, 1:L[h]])
  }
  S1 <- array(dim = c(n12, p))
  S2 <- array(dim = c(n12, p))

  for(i in 1:n12){
    for (h in 1:p) {
      S1[i, h] <- rbinom(n = 1, size = 1, prob = beta[h])
      if(S1[i, h] == 0){
        X1[i, h] <- X12[i, h]
      }else{
        X1[i, h] <- sample.int(n = L[h], size = 1, prob = theta[h, 1:L[h]])
      }

      S2[i, h] <- rbinom(n = 1, size = 1, prob = beta[h])
      if(S2[i, h] == 0){
        X2[i, h] <- X12[i, h]
      }else{
        X2[i, h] <- sample.int(n = L[h], size = 1, prob = theta[h, 1:L[h]])
      }
    }
  }

  Delta <- Matrix::sparseMatrix(c(1:n12), c(1:n12), x = 1, dims = c(n1, n2))
  Delta <- as(Delta, "TsparseMatrix")

  return(list("X1" = X1,
              "X2" = X2,
              "beta" = beta,
              "theta" = theta,
              "L" = L,
              "Delta" = Delta
              ))

}
