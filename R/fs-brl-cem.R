
#' @rdname brl_cem
#' @export
fs_brl_cem <- function(X1,
                       X2,
                       a_prop = 2,
                       b_prop = 2,
                       max_iter = 50,
                       eps = 1e-4,
                       complete_matching = FALSE,
                       candidate_pairs = NULL,
                       comp_data = NULL,
                       solver = "clp") {


  sol <- resolve_roi_solver(solver)
  solver <- sol$solver
  control_list <- sol$control

  if (is.null(comp_data)) {
    prep <- prep_brl_inputs(X1 = X1,
                            X2 = X2,
                            candidate_pairs = candidate_pairs)
    X1 <- prep$X1
    X2 <- prep$X2
    n1 <- prep$n1
    n2 <- prep$n2
    p <- prep$p
    candidate_pairs <- prep$candidate_pairs
    switched_order <- prep$switched_order

    out_recode <- recode_columns(X1, X2)
    X1 <- out_recode$X1
    X2 <- out_recode$X2
    Gamma <- compare_binary(X1, X2, how = "binary", candidate_pairs = candidate_pairs)
    out_hash <- hash_Gamma(Gamma = Gamma)
  }else{
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    p <- ncol(X1)
    switched_order <- FALSE
    Gamma <- comp_data$Gamma
    out_hash <- comp_data$out_hash

    if (!is.null(candidate_pairs)) {
      lin_idx_all <- (candidate_pairs[, 2] - 1L) * n1 + candidate_pairs[, 1]
      if (!all(comp_data$out_hash$lin_idx_nz %in% lin_idx_all)) {
        stop("candidate_pairs do not match the comp_data that was supplied.")
      }
    }
  }

  if(isFALSE(complete_matching)){
    dir_lp_constraints <- rep("<=", n1 + n2)
  }else{
    dir_lp_constraints <- c(rep("<=", n1), rep("==", n2))
  }

  if(is.null(candidate_pairs)){
    n_candidate_pairs <- n1*n2
  }else{
    n_candidate_pairs <- nrow(candidate_pairs)

    lin_idx_all <- (candidate_pairs[, 2] - 1L) * n1 + candidate_pairs[, 1]
    match_idx <- match(lin_idx_all, out_hash$lin_idx_nz, nomatch = 0L)
    pairs_code_all <- integer(length(lin_idx_all))
    pairs_code_all[match_idx != 0L] <- out_hash$pairs_code_nz[match_idx[match_idx != 0L]]
    pairs_code_all[match_idx == 0L] <- 1L
  }

  # parameter initialization
  L <- rep(2, p)
  freq_Gamma <- vector(mode = "list", length = p)
  for(h in 1:p){
    freq_Gamma[[h]] <- c(n_candidate_pairs - length(Gamma[[h]]@x),
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

    if (!is.null(candidate_pairs)) {
      if (isFALSE(complete_matching)) {
        which_codes <- which(cost_unique > 0)
        keep <- pairs_code_all %in% which_codes
        idx_long_filtered <- candidate_pairs[keep, , drop = FALSE]
        cost_vec <- cost_unique[pairs_code_all[keep]]
      }else{
        cost_vec <- cost_unique[pairs_code_all]
        idx_long_filtered <- candidate_pairs
      }
    } else {
      which_codes <- which(cost_unique > 0)
      keep_nz <- out_hash$pairs_code_nz %in% which_codes
      which_pairs <- out_hash$lin_idx_nz[keep_nz]
      cost_vec <- cost_unique[out_hash$pairs_code_nz[keep_nz]]
      i_idx <- ((which_pairs - 1L) %% n1) + 1L
      j_idx <- ((which_pairs - 1L) %/% n1) + 1L
      idx_long_filtered <- cbind(i_idx, j_idx)
    }

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
                                                        dir = dir_lp_constraints,
                                                        rhs = rep(1, n1 + n2)),
                        maximum = TRUE)
      out_roi <- ROI::ROI_solve(lp_obj, solver = solver, control = control_list)
      if (!isTRUE(out_roi$status$code == 0) || length(out_roi$solution) == 0) {
        stop("LP solver failed or returned no solution.")
      }
      if (isTRUE(complete_matching) && sum(out_roi$solution == 1) < n2) {
        stop("Infeasible complete matching with the supplied candidate graph.")
      }
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

  id_X2 <- get_id_X2(n1, n2, coreferent_pairs)

  out <- list("Delta" = Delta,
              "m" = m,
              "u" = u,
              "prop" = prop,
              "coreferent_pairs" = coreferent_pairs,
              "id_X2" = id_X2,
              "cloglik" = cloglik[1:(em_iter-1)],
              "switched_order" = switched_order)

  return(out)

}
