
# NOT USED: setting the Lagrangian equal to 0

# ftheta_lagrange <- function(lambda, beta, n12, nx){
#   a <- lambda * beta * (2 - beta)
#   b <- beta * (2 - beta) * nx + lambda * (1 - beta) ^ 2
#   c <- (1 - beta) ^ 2 * (nx - n12)
#   return(-(-b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a))
# }

# sum_theta_lagrange <- function(lambda, beta, n12, nx) {
#   sum(ftheta_lagrange(lambda, beta, n12, nx)) - 1
# }


#################################################################
##                       Score functions                       ##
#################################################################

loglik_theta_beta_h <- function(theta_h, beta_h, n12_h, nx_h, a_beta, b_beta){
  sum(n12_h[-length(n12_h)] * log(beta_h * (2 - beta_h) + (1 - beta_h) ^ 2 / theta_h)) + n12_h[length(n12_h)]*log(beta_h*(2 - beta_h)) + sum(nx_h * log(theta_h)) +
    (a_beta - 1)*log(beta_h) + (b_beta - 1)*log(1 - beta_h)
}

score_fun_beta <- function(beta, theta, n12_h, a_beta, b_beta){
  sum((2 * (1 - beta) * (1 - 1 / theta[theta > 0])) / (beta * (2 - beta) + (1 - beta) ^ 2 / theta[theta > 0]) * n12_h[-length(n12_h)][theta > 0])+
    (2 * (1 - beta)) / (beta * (2 - beta)) * n12_h[length(n12_h)] +
    (a_beta - 1)/beta - (b_beta - 1)/(1 - beta)
}



##################################################################
##                       Numerical M-step                       ##
##################################################################

softmax <- function(w) {
  exp_w <- exp(w - max(w))
  return(exp_w / sum(exp_w))
}

log_softmax <- function(w) {
  w - max(w) - log(sum(exp(w - max(w))))
}

loglik_h <- function(params, n12_h, nx_h, a_beta, b_beta){
  q <- params[1]
  w <- params[-1]
  -(
    sum(n12_h[-length(n12_h)] * log(plogis(q) * (2 - plogis(q)) + (1 - plogis(q)) ^ 2 / softmax(w)))
    + n12_h[length(n12_h)]*log(plogis(q)*(2 - plogis(q))) + sum(nx_h * log(softmax(w))) +
      dbeta(plogis(q), a_beta, b_beta, log = TRUE)
  )
}

m_step_h_optim <- function(beta, theta, n12_h, nx_h, a_beta, b_beta){
  L_h <- length(nx_h)
  out_optim_h <- optim(par = c(qlogis(beta), log(theta)),
                       fn = loglik_h,
                       n12_h = n12_h,
                       nx_h = nx_h,
                       a_beta = a_beta,
                       b_beta = b_beta)
  beta <- plogis(out_optim_h$par[1])
  theta <- softmax(out_optim_h$par[-1])
  return(list("beta" = beta,
              "theta" = theta))
}


#' @rdname brl_cem
#' @export
graph_brl_cem <- function(X1,
                          X2,
                          a_prop = 2,
                          b_prop = 2,
                          a_beta = 1.5,
                          b_beta = 5,
                          theta_method = c("empirical", "uniform", "estimated"),
                          max_iter = 50,
                          eps = 1e-4,
                          candidate_pairs = NULL,
                          comp_data = NULL,
                          solver = "clp") {

  sol <- resolve_roi_solver(solver)
  solver <- sol$solver
  control_list <- sol$control

  theta_method <- match.arg(theta_method)

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

  L <- apply(rbind(X1, X2), 2, \(x) length(unique(x)))

  if(is.null(comp_data)){
    out_recode <- recode_columns(X1, X2)
    X1 <- out_recode$X1
    X2 <- out_recode$X2
    is_binary_comp <- theta_method == "uniform"
    Gamma <- compare_binary(X1, X2, how = ifelse(is_binary_comp, "binary", "value"), candidate_pairs = candidate_pairs)
    out_hash <- hash_Gamma(Gamma = Gamma, L = ifelse(is_binary_comp, NULL, L))
  }else{
    Gamma <- comp_data$Gamma
    out_hash <- comp_data$out_hash
  }

  theta <- array(dim = c(p, max(L)))
  nx <- array(dim = dim(theta))
  for(h in 1:p){
    nx[h, 1:L[h]] <- tabulate(bin = c(X1[, h], X2[, h]), nbins = L[h])
    if(theta_method == "uniform"){
      theta[h, 1:L[h]] <- rep(1/L[h], L[h])
    }else{
      theta[h, 1:L[h]] <- (nx[h, 1:L[h]] + 1)/(n1 + n2 + L[h])
    }
  }

  # initialization
  beta <- runif(n = p, min = 0.01, max = 0.4)
  prop <- rbeta(n = 1, shape1 = 7, shape2 = 2)

  if(is.null(candidate_pairs)){
    # candidate_pairs <- cbind(rep(1:n1, times = n2), rep(1:n2, each = n1))
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

    if(theta_method == "uniform"){
      cost_unique <- numeric(nrow(out_hash$Gamma_unique))
      for(i in 1:length(cost_unique)){
        cost_unique[i] <- sum(log(beta * (2 - beta) + (1 - beta) ^ 2 * L * (out_hash$Gamma_unique[i, ] - 1))) + log(prop) - log(n1) - log(1 - prop)
      }
    }else{
      theta_extra <- cbind(Inf, theta)
      cost_unique <- numeric(nrow(out_hash$Gamma_unique))
      for(i in 1:length(cost_unique)){
        cost_unique[i] <- sum(log(beta * (2 - beta) + (1 - beta) ^ 2 / theta_extra[cbind(1:p, out_hash$Gamma_unique[i, ])])) + log(prop) - log(n1) - log(1 - prop)
      }
    }
    which_codes <- which(cost_unique > 0)

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
      out_roi <- ROI::ROI_solve(lp_obj, solver = solver, control = control_list)
      coreferent_pairs <- idx_long_filtered[out_roi$solution == 1, , drop = FALSE]

      cloglik[em_iter] <- out_roi$objval
    }

    cloglik[em_iter] <- cloglik[em_iter] + n2 * log(1 - prop) +
      dbeta(x = prop, shape1 = a_prop, shape2 = b_prop, log = TRUE) +
      sum(dbeta(x = beta, shape1 = a_beta, shape2 = b_beta, log = TRUE))
    for(h in 1:p){
      cloglik[em_iter] <- cloglik[em_iter] + sum(log(theta[h, c(X1[, h], X2[, h])]))
    }

    cat(em_iter, "/ ")

    if(em_iter > 1){
      cll_increase <- cloglik[em_iter] - cloglik[em_iter-1]
    }
    em_iter <- em_iter + 1

    for(h in 1:p){
      x12_h <- ifelse(X1[coreferent_pairs[, 1], h] == X2[coreferent_pairs[, 2], h], X1[coreferent_pairs[, 1], h], L[h] + 1)
      if(length(x12_h) == 0){
        n12_h <- rep(0, L[h] + 1)
      }else{
        n12_h <- tabulate(x12_h, L[h] + 1)
      }
      if(theta_method != "estimated"){

        # beta[h] <- optimize(
        #   f = function(x) -loglik_theta_beta_h(theta[h, 1:L[h]], x, n12_h, nx[h, 1:L[h]], a_beta, b_beta),
        #   interval = c(1e-6, 1 - 1e-6)
        # )$minimum

        beta[h] <- tryCatch(expr = uniroot(f = score_fun_beta,
                                           interval = c(1e-6, 1-1e-6),
                                           theta = theta[h, 1:L[h]],
                                           n12_h = n12_h,
                                           a_beta = a_beta,
                                           b_beta = b_beta)$root,
                            error = function(cond) {
                              ifelse(loglik_theta_beta_h(theta[h, 1:L[h]], 1e-6, n12_h, nx[h, 1:L[h]], a_beta, b_beta) > loglik_theta_beta_h(theta[h, 1:L[h]], 1-1e-6, n12_h, nx[h, 1:L[h]], a_beta, b_beta),
                                     1e-6,
                                     1 - 1e-6)}
                            )
      }else{
        out_m_step_h <- m_step_h_optim(beta = beta[h], theta = theta[h, 1:L[h]], n12_h = n12_h, nx_h = nx[h, 1:L[h]],
                                       a_beta = a_beta,
                                       b_beta = b_beta)
        theta[h, 1:L[h]] <- out_m_step_h$theta
        beta[h] <- out_m_step_h$beta
      }
    }

    n12 <- nrow(coreferent_pairs)
    prop <- (n12 + a_prop - 1) / (n2 + a_prop + b_prop - 2)

  }

  Delta <- Matrix::sparseMatrix(i = coreferent_pairs[, 1], coreferent_pairs[, 2], x = 1, dims = c(n1, n2))
  Delta <- as(Delta, "TsparseMatrix")

  out <- list("Delta" = Delta,
              "beta" = beta,
              "theta" = theta,
              "prop" = prop,
              "coreferent_pairs" = coreferent_pairs,
              "cloglik" = cloglik[1:(em_iter-1)],
              "theta_method" = theta_method,
              "switched_order" = switched_order)

  return(out)
}

