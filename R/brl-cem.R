#' CEM estimation for bipartite record linkage models
#'
#' Fit bipartite record linkage models via a Classification EM (CEM)
#' algorithm. The matching problem is solved using
#' \pkg{ROI} with the \pkg{ROI.plugin.clp} backend.
#'
#' Two models are supported:
#' \itemize{
#'   \item \code{fs_brl_cem()}: the Fellegi–Sunter model (two-class latent model fitted to pairwise comparison data),
#'   \item \code{graph_brl_cem()}: the graphical record linkage.
#' }
#' Both models use the \emph{fast beta linkage} (fabl) prior
#' for the matching structure.
#'
#' The function \code{brl_cem()} is a wrapper that runs either model multiple
#' times and returns the fit with the highest final complete data log-likelihood.
#'
#' @param X1,X2 Matrices of the two categorical datasets to be linked.
#' Rows correspond to records and columns to linking variables.
#' Both must have the same variables in the same order.
#'
#' @param model Character scalar used by \code{brl_cem()}:
#'   \code{"fs"} selects the Fellegi–Sunter model,
#'   \code{"graph"} selects the graphical record linkage model.
#'
#' @param reps Integer. Number of repeated CEM fits in \code{brl_cem()}.
#'
#' @param candidate_pairs Optional integer matrix with two columns specifying
#' candidate record pairs to be considered for matching. Each row corresponds
#' to a pair \code{(i, j)}, where \code{i} indexes rows of \code{X1} and
#' \code{j} indexes rows of \code{X2} (both 1-based).
#' If \code{NULL} (default), all possible record pairs are considered.
#'
#' @param a_prop,b_prop Positive numeric scalars specifying Beta prior parameters
#'   for the overlap proportion.
#'
#' @param a_beta,b_beta Positive numeric scalars specifying Beta prior parameters
#'   for the distortion probabilities.
#'
#' @param theta_method Character scalar controlling how category probabilities
#' are handled in the graphical model. One of \code{"empirical"},
#' \code{"uniform"}, or \code{"estimated"}.
#'
#' @param max_iter Integer. Maximum number of CEM iterations.
#'
#' @param eps Numeric. Convergence tolerance on the increase in complete
#'   log-likelihood.
#'
#' @param ... Additional arguments passed to the selected fitting function when
#' using \code{brl_cem()}.
#'
#' @param comp_data Optional list of precomputed comparison data.
#'
#' @param solver Character string specifying the LP solver to be used by
#'   \pkg{ROI}. Currently supported options are \code{"clp"} and
#'   \code{"lpsolve"}. The default is \code{"clp"}; if the corresponding
#'   plugin package \pkg{ROI.plugin.clp} is not available, a warning is
#'   issued and the solver falls back to \code{"lpsolve"} (if available).
#'   If \code{"lpsolve"} is requested explicitly, the package
#'   \pkg{ROI.plugin.lpsolve} must be installed.
#'
#' @return
#' Each routine returns a list containing:
#' \describe{
#'   \item{Delta}{A sparse \code{TsparseMatrix} encoding the estimated bipartite
#'     matching between records in \code{X1} and \code{X2}.}
#'   \item{cloglik}{Numeric vector of complete data log-likelihood values across CEM
#'     iterations.}
#'   \item{prop}{Estimated overlap proportion of linked records.}
#'   \item{coreferent_pairs}{Matrix with indices of the predicted coreferent pairs.}
#'   \item{switched_order}{Logical; \code{TRUE} if the input datasets were swapped
#'     internally because \code{nrow(X1) < nrow(X2)}.}
#' }
#'
#' Model-specific outputs:
#' \itemize{
#'   \item \code{fs_brl_cem()}: returns variable-specific agreement probabilities for
#'   coreferent and non coreferent pairs \code{m} and \code{u}.
#'   \item \code{graph_brl_cem()}: returns distortion probabilities \code{beta}
#'     and category probabilities \code{theta}, and the chosen \code{theta_method}.
#'   \item \code{brl_cem()}: additionally returns the number of runs of the CEM algorithm
#'   (\code{reps}) and the number of failed runs (\code{error_reps}).
#' }
#'
#' @references
#' Kundinger, B., Reiter, J. P., and Steorts, R. C. (2025).
#' \emph{Efficient and Scalable Bipartite Matching with Fast Beta Linkage (fabl)}.
#' Bayesian Analysis, 20(3), 949--972.
#' \doi{10.1214/24-BA1427}
#'
#' Steorts, R. C., Hall, R., and Fienberg, S. E. (2016).
#' \emph{A Bayesian Approach to Graphical Record Linkage and Deduplication}.
#' Journal of the American Statistical Association, 111(516), 1660--1672.
#' \doi{10.1080/01621459.2015.1105807}
#'
#' @examples
#' \dontrun{
#' # Fellegi–Sunter model
#' fit_fs <- fs_brl_cem(X1, X2, max_iter = 20)
#'
#' # Graphical record linkage model
#' fit_graph <- graph_brl_cem(X1, X2, max_iter = 20)
#'
#' # Repeated runs, best solution retained
#' fit_reps <- brl_cem(X1, X2, model = "graph", reps = 5, max_iter = 20)
#' }
#'
#' @name brl_cem
NULL

#' @rdname brl_cem
#' @export
brl_cem <- function(X1,
                    X2,
                    model = c("fs", "graph"),
                    reps = 5,
                    candidate_pairs = NULL,
                    solver = "clp",
                    ...) {

  sol <- resolve_roi_solver(solver)
  solver <- sol$solver

  model <- match.arg(model)

  foo <- switch(model,
                "graph" = graph_brl_cem,
                "fs" = fs_brl_cem)

  extra_args <- list(...)
  if (is.null(extra_args$theta_method) & model == "graph") {
    extra_args$theta_method <- "empirical"
  }

  out_recode <- recode_columns(X1, X2)
  X1 <- out_recode$X1
  X2 <- out_recode$X2

  if(model == "graph" && extra_args$theta_method != "uniform"){
    L <- apply(rbind(X1, X2), 2, \(x) length(unique(x)))
  }else{
    L <- NULL
  }
  how <- switch(model,
                "graph" = ifelse(extra_args$theta_method == "uniform", "binary", "value"),
                "fs"    = "binary")

  comp_data <- list()
  comp_data$Gamma <- compare_binary(X1, X2, how, candidate_pairs)
  comp_data$out_hash <- hash_Gamma(Gamma = comp_data$Gamma, L = L)

  out <- vector(mode = "list", length = reps)

  for(i in seq_len(reps)){
    cat("\nModel:", i, "\n\t")

    out[[i]] <- tryCatch(
      expr = do.call(foo, c(list(X1 = X1, X2 = X2, comp_data = comp_data, candidate_pairs = candidate_pairs, solver = solver), extra_args)),
      error = function(e) {
        message("Error caught: ", e$message)
        return(NA)
      }
    )
  }

  is_error <- vapply(out, function(x) identical(x, NA), logical(1))
  n_error  <- sum(is_error)

  if(n_error == reps){
    stop("All repetitions failed.")
  }
  out <- out[!is_error]
  best_rep <- which.max(sapply(out, \(x) tail(x$cloglik, 1)))

  out <- out[[best_rep]]
  out$reps <- reps
  out$error_reps <- n_error

  return(out)
}
