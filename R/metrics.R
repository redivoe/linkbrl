
#' Binary classification metrics
#'
#' Compute precision, recall, and F1 score for binary predictions.
#'
#' The function supports both dense binary objects (vectors, matrices)
#' and sparse matrices of class \code{dgTMatrix}. In the sparse case,
#' non-zero entries are interpreted as positive class predictions.
#'
#' @param Delta A binary vector, matrix, or sparse matrix of class
#'   \code{dgTMatrix} representing the true labels.
#'
#' @param Delta_hat A binary vector, matrix, or sparse matrix of class
#'   \code{dgTMatrix} representing the predicted labels.
#'
#' @return A named numeric vector with elements:
#'   \describe{
#'     \item{recall}{True positive rate.}
#'     \item{precision}{Positive predictive value.}
#'     \item{F1}{F1 score (harmonic mean of precision and recall).}
#'   }
#'
#' @examples
#' # Dense case
#' y <- c(1, 0, 1, 1, 0)
#' y_hat <- c(1, 0, 0, 1, 1)
#' binary_metrics(y, y_hat)
#'
#' @export
binary_metrics <- function(Delta, Delta_hat){

  Delta_has_dim <- !is.null(dim(Delta))
  Delta_hat_has_dim <- !is.null(dim(Delta_hat))
  stopifnot("Delta and Delta_hat must both be vectors or both be matrices" =
              Delta_has_dim == Delta_hat_has_dim)

  if(Delta_has_dim){
    stopifnot("The two matrices have different dimensions." =
                identical(dim(Delta), dim(Delta_hat)))
  } else {
    stopifnot("The two vectors are of different length." =
                length(Delta) == length(Delta_hat))
  }

  if(inherits(Delta, "dgTMatrix") & inherits(Delta_hat, "dgTMatrix")){
    P <- length(Delta@x)
    PP <- length(Delta_hat@x)
    TP <- length(intersect(x = paste(Delta@i, Delta@j),
                           y = paste(Delta_hat@i, Delta_hat@j)))
  }else{
    y <- as.vector(Delta)
    y_hat <- as.vector(Delta_hat)
    P <- sum(y == 1)
    PP <- sum(y_hat == 1)
    TP <- sum(y == 1 & y_hat == 1)
  }

  precision <- if (PP == 0) {
    0
  } else{
    TP / PP
  }

  recall <- if (P == 0) {
    0
  } else{
    TP / P
  }

  F1 <- if (precision + recall == 0) {
    0
  } else{
    2 * precision * recall / (precision + recall)
  }

  return(c("recall" = recall, "precision" = precision, "F1" = F1))
}
