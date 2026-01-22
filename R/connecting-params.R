muparams_graph_er <- function(beta, theta){
  p <- length(beta)
  u <- rowSums(theta^2, na.rm = TRUE)
  m <- (1 - beta)^2 +
    2 * (1 - beta) * beta * rowSums(theta^2, na.rm = TRUE) +
    beta^2 * rowSums(theta^2, na.rm = TRUE)

  return(list("m" = m,
              "u" = u))
}
