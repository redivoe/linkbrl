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
