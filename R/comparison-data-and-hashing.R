
compare_binary <- function(X1, X2, how, candidate_pairs = NULL){

  compare_fun <- switch(how,
                        "binary" = `==`,
                        "value"  = \(x, y) ifelse(x == y, x, 0),
                        stop("Unknown 'how' in compare_binary()."))

  p  <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  if (is.null(candidate_pairs)) {
    i <- rep(seq_len(n1), times = n2)
    j <- rep(seq_len(n2), each = n1)
  } else {
    if (!is.matrix(candidate_pairs) || ncol(candidate_pairs) != 2) {
      stop("candidate_pairs must be an (n_candidate_pairs x 2) integer matrix.")
    }
    i <- as.integer(candidate_pairs[, 1])
    j <- as.integer(candidate_pairs[, 2])
  }

  Gamma <- vector("list", p)
  cat("\nComparing datasets: ")

  for (h in seq_len(p)) {
    x <- compare_fun(X1[i, h], X2[j, h])

    if (how == "binary") {
      nz     <- x
      x_vals <- rep.int(1, sum(nz))
    } else {  # how == "value"
      nz     <- x != 0
      x_vals <- as.double(x[nz])
    }

    if (any(nz)) {
      Gamma[[h]] <- Matrix::sparseMatrix(i = i[nz],
                                         j = j[nz],
                                         x = x_vals,
                                         dims = c(n1, n2))
    } else {
      Gamma[[h]] <- Matrix::sparseMatrix(i = integer(0),
                                         j = integer(0),
                                         x = numeric(0),
                                         dims = c(n1, n2))
    }
    cat(".")
  }
  cat("/ ")
  attr(Gamma, "how") <- how
  return(Gamma)
}


hash_Gamma <- function(Gamma, L = NULL){
  p <- length(Gamma)

  is_binary <- attr(Gamma, "how") == "binary"

  if(is_binary){
    if(!is.null(L)){
      stop("L is provided to a binary Gamma.")
    }
    L <- rep(2, p)
  }else if(attr(Gamma, "how") == "value"){
    if(is.null(L)){
      stop("L must be provided to a value Gamma.")
    }
  }else{
    stop("Gamma provided to hash_Gamma of unknown type.")
  }

  order_sparse <- order(sapply(Gamma, \(x) length(x@x)))
  Gamma <- Gamma[order_sparse]

  if(is_binary){
    L_cumsum <- cumsum(L[order_sparse] - 1)

    hash <- Gamma[[1]]
    hash@x <- (2^(L_cumsum[1]))[hash@x]
    for(h in 2:p){
      temp <- Gamma[[h]]
      temp@x <- (2^(L_cumsum[h]))[temp@x]
      hash <- hash + temp
    }

  }else{
    hash <- Gamma[[1]]
    for (h in 2:p) {
      hash <- hash * 31 + Gamma[[h]]
    }
  }

  hash <- Matrix::Matrix(hash, sparse = TRUE)
  hash <- methods::as(hash, "TsparseMatrix")

  n1 <- hash@Dim[1]
  ij <- cbind(hash@i + 1L, hash@j + 1L)
  lin_idx <- (ij[,2] - 1L) * n1 + ij[,1]

  not_dupe <- !duplicated(hash@x)
  unique_indices <- ij[not_dupe, , drop = FALSE]

  hash@x <- as.double(match(hash@x, hash@x[not_dupe]))
  pairs_code_nz <- as.integer(hash@x) + 1L

  Gamma <- Gamma[order(order_sparse)]
  Gamma_unique <- matrix(nrow = nrow(unique_indices), ncol = p)
  for(h in 1:p){
    Gamma_unique[, h] <- Gamma[[h]][unique_indices]
  }
  Gamma_unique <- rbind(rep(0, p),
                        Gamma_unique) + 1

  return(list(
    "pairs_code_nz" = pairs_code_nz,
    "lin_idx_nz" = lin_idx,
    "Gamma_unique" = Gamma_unique
  ))
}

