recode_columns <- function(X1, X2) {
  p <- ncol(X1)
  X1_new <- X1
  X2_new <- X2

  for (h in seq_len(p)) {
    vals <- sort(unique(c(X1[, h], X2[, h])))
    X1_new[, h] <- match(X1[, h], vals)
    X2_new[, h] <- match(X2[, h], vals)
  }

  return(list(X1 = X1_new, X2 = X2_new))
}


resolve_roi_solver <- function(solver = c("clp", "lpsolve")) {
  solver <- match.arg(solver)

  has_clp <- requireNamespace("ROI.plugin.clp", quietly = TRUE)
  has_lpsolve <- requireNamespace("ROI.plugin.lpsolve", quietly = TRUE)

  if (solver == "clp" && !has_clp) {
    warning("Package {ROI.plugin.clp} not available. Falling back to solver = 'lpsolve'.")
    solver <- "lpsolve"
  }

  if (solver == "lpsolve" && !has_lpsolve) {
    stop("Package {ROI.plugin.lpsolve} is required when solver = 'lpsolve'.", call. = FALSE)
  }

  control <- if (solver == "clp") list(amount = 0) else list()

  return(list(solver = solver, control = control))
}
