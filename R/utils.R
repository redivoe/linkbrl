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

prep_brl_inputs <- function(X1, X2, candidate_pairs = NULL) {

  # basic checks
  if (ncol(X1) != ncol(X2)) {
    stop("X1 and X2 must have the same number of columns (referring to the same variables in the same order).")
  }
  if (anyNA(rbind(X1, X2))) {
    stop("No missing values must be present in either X1 or X2.")
  }

  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p  <- ncol(X1)

  # validate candidate_pairs in original orientation
  if (!is.null(candidate_pairs)) {
    if (!is.matrix(candidate_pairs) || ncol(candidate_pairs) != 2L) {
      stop("candidate_pairs must be an (n x 2) matrix.")
    }
    storage.mode(candidate_pairs) <- "integer"
    if (any(candidate_pairs < 1L)) {
      stop("candidate_pairs contains indices < 1.")
    }
    if (any(candidate_pairs[, 1] > n1) || any(candidate_pairs[, 2] > n2)) {
      stop("candidate_pairs contains indices out of bounds for X1/X2.")
    }
  }

  # swap so that n2 <= n1
  switched_order <- FALSE
  if (n2 > n1) {
    switched_order <- TRUE

    tmp <- X2
    X2 <- X1
    X1 <- tmp

    n1 <- nrow(X1)
    n2 <- nrow(X2)

    if (!is.null(candidate_pairs)) {
      candidate_pairs <- candidate_pairs[, c(2, 1), drop = FALSE]
    }
  }

  return(
    list(
      X1 = X1,
      X2 = X2,
      n1 = n1,
      n2 = n2,
      p = p,
      candidate_pairs = candidate_pairs,
      switched_order = switched_order
    )
  )
}

get_id_X2 <- function(n1, n2, coreferent_pairs){
  id_X2 <- integer(n2)
  if (nrow(coreferent_pairs) > 0) {
    id_X2[coreferent_pairs[, 2]] <- coreferent_pairs[, 1]
  }
  which_not_matched_X2 <- setdiff(seq_len(n2), coreferent_pairs[, 2])
  if (length(which_not_matched_X2) > 0) {
    id_X2[which_not_matched_X2] <- n1 + seq_len(length(which_not_matched_X2))
  }
  return(id_X2)
}

