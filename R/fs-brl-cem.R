
#' @rdname brl_cem
#' @export
fs_brl_cem <- function(X1,
                       X2,
                       a_prop = 2,
                       b_prop = 2,
                       max_iter = 50,
                       eps = 1e-4,
                       candidate_pairs = NULL,
                       comp_data = NULL) {

  if (!requireNamespace("ROI.plugin.clp", quietly = TRUE)) {
    stop("Package ROI.plugin.clp is needed.")
  }

  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)
  switched_order <- FALSE

  if(n2 > n1){
    switched_order <- TRUE
    temp <- X2
    X2 <- X1
    X1 <- temp
    n1 <- nrow(X1)
    n2 <- nrow(X2)
  }

  if(is.null(comp_data)){
    out_recode <- recode_columns(X1, X2)
    X1 <- out_recode$X1
    X2 <- out_recode$X2
    Gamma <- compare_binary(X1, X2, how = "binary", candidate_pairs = candidate_pairs)
    out_hash <- hash_Gamma(Gamma = Gamma)
  }else{
    Gamma <- comp_data$Gamma
    out_hash <- comp_data$out_hash
  }

  L <- rep(2, p)
  freq_Gamma <- vector(mode = "list", length = p)
  for(h in 1:p){
    freq_Gamma[[h]] <- c(n1*n2 - length(Gamma[[h]]@x),
                         tabulate(Gamma[[h]]@x, L[h] - 1))
  }

  alpha_m <- lapply(L, \(L_h) seq(1.5, 3, len = L_h))
  alpha_u <- lapply(L, \(L_h) seq(3, 1.5, len = L_h))

  m <- array(dim = c(p, max(L)))
  u <- array(dim = c(p, max(L)))
  if(all(L == 2)){
    m[, 2] <- runif(n = p, min = 0.6, max = 0.99)
    u[, 2] <- runif(n = p, min = 0.01, max = 0.4)
    m[, 1] <- 1 - m[, 2]
    u[, 1] <- 1 - u[, 2]
  }else{
    for(h in 1:p){
      m[h, 1:L[h]] <- MCMCpack::rdirichlet(n = 1, alpha = alpha_m[[h]])
      u[h, 1:L[h]] <- MCMCpack::rdirichlet(n = 1, alpha = alpha_u[[h]])
    }
  }
  prop <- rbeta(n = 1, shape1 = 7, shape2 = 2)

  if(is.null(candidate_pairs)){
    n_candidate_pairs <- n1*n2
  }else{
    n_candidate_pairs <- nrow(candidate_pairs)
  }
  coreferent_pairs <- matrix(nrow = 0, ncol = 2)


  cloglik <- rep(0, max_iter)
  cll_increase <- 1

  em_iter <- 1
  cat("Iteration: ")

  while((em_iter <= max_iter) & (cll_increase > eps)){

    cost_unique <- numeric(nrow(out_hash$Gamma_unique))
    for(i in 1:length(cost_unique)){
      idx <- cbind(1:p, out_hash$Gamma_unique[i, ])
      cost_unique[i] <- sum(log(m[idx])) - sum(log(u[idx])) + log(prop) - log(n1) - log(1 - prop)
    }
    which_codes <- which(cost_unique > 0)

    # which_pairs <- which(out_hash$pairs_code %in% which_codes)
    # cost_vec <- cost_unique[out_hash$pairs_code[which_pairs]]

    keep_nz <- out_hash$pairs_code_nz %in% which_codes
    which_pairs <- out_hash$lin_idx_nz[keep_nz]
    cost_vec <- cost_unique[out_hash$pairs_code_nz[keep_nz]]
    i_idx <- ((which_pairs - 1L) %% n1) + 1L
    j_idx <- ((which_pairs - 1L) %/% n1) + 1L
    idx_long_filtered <- cbind(i_idx, j_idx)

    # idx_long_filtered <- candidate_pairs[which_pairs, , drop = FALSE]
    if(nrow(idx_long_filtered) == 0){
      warning("No pairs matchable (with positive cost)")
    }else{
      constr_stm <- slam::simple_triplet_matrix(i = c(idx_long_filtered[, 1], idx_long_filtered[, 2] + n1),
                                                j = rep(1:nrow(idx_long_filtered), 2),
                                                v = rep(1, 2 * nrow(idx_long_filtered)),
                                                nrow = n1 + n2,
                                                ncol = nrow(idx_long_filtered))

      lp_obj <- ROI::OP(objective = cost_vec,
                        constraints = ROI::L_constraint(L = constr_stm,
                                                        dir = rep('<=', n1 + n2),
                                                        rhs = rep(1, n1 + n2)),
                        maximum = TRUE)
      out_roi <- ROI::ROI_solve(lp_obj, solver = "clp", control = list(amount = 0))
      coreferent_pairs <- idx_long_filtered[out_roi$solution == 1, , drop = FALSE]

      cloglik[em_iter] <- out_roi$objval
    }

    cloglik[em_iter] <- cloglik[em_iter] + n2 * log(1 - prop) +
      dbeta(x = prop, shape1 = a_prop, shape2 = b_prop, log = TRUE)
    for(h in 1:p){
      cloglik[em_iter] <- cloglik[em_iter] +
        sum(log(u[h, ]) * freq_Gamma[[h]]) +
        log(MCMCpack::ddirichlet(x = m[h, 1:L[h]], alpha = alpha_m[[h]])) +
        log(MCMCpack::ddirichlet(x = u[h, 1:L[h]], alpha = alpha_u[[h]]))
    }

    cat(em_iter, "/ ")

    if(em_iter > 1){
      cll_increase <- cloglik[em_iter] - cloglik[em_iter-1]
    }
    em_iter <- em_iter + 1

    n12 <- nrow(coreferent_pairs)
    # M-step
    for(h in 1:p){
      if(n12 == 0){
        freq_12h <- rep(0, L[h])
      }else{
        freq_12h <- tabulate(Gamma[[h]][coreferent_pairs] + 1, L[h])
      }
      m[h, 1:L[h]] <- (freq_12h + alpha_m[[h]] - 1)/(n12 + sum(alpha_m[[h]]) - L[h])
      u[h, 1:L[h]] <- (freq_Gamma[[h]] - freq_12h + alpha_u[[h]] - 1)/(n_candidate_pairs - n12 + sum(alpha_u[[h]]) - L[h])
      u[h, 1:L[h]] <- u[h, 1:L[h]] / sum(u[h, 1:L[h]])
    }
    prop <- (n12 + a_prop - 1) / (n2 + a_prop + b_prop - 2)

  }

  Delta <- Matrix::sparseMatrix(i = coreferent_pairs[, 1], coreferent_pairs[, 2], x = 1, dims = c(n1, n2))
  Delta <- as(Delta, "TsparseMatrix")

  out <- list("Delta" = Delta,
              "m" = m,
              "u" = u,
              "prop" = prop,
              "coreferent_pairs" = coreferent_pairs,
              "cloglik" = cloglik[1:(em_iter-1)],
              "switched_order" = switched_order)

  return(out)

}
