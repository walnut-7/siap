# Required packages:
# - Matrix
## If not included, 'Error in t.default(M) : argument is not a matrix.'
siap <- function(
    X, # If there are NA's in X, X is considered partially observed
    # --------- step 1 ----------
    rank.max,
    step1 = T,
    step1.prd = T,
    step1.Z = NULL, # optional generic mean design for phase 1; ncol(X) x q
    use.ecm = T,
    incl.cov = NULL,
    lambda = 0,
    lambda_mu1 = 0,
    # --------- step 2 ----------
    fit1 = NULL, # output of step 1
    step2.rank = 10,
    step2 = T,
    step2.ar = T,
    step2.prd = T,
    step2.Z = NULL, # optional generic mean design for phase 2; ncol(X) x q
    step2.b.update = c("iterative", "whole", "auto"), # strategy for updating b in phase 2
    lambda1 = 0,
    lambda2 = 0,
    lambda_mu2 = 0,
    alpha = 0,
    time.lag = 2,
    # ---------------- solving strategy ---------------
    beta = NULL, # AR parameters. A p x r matrix, indicating the diagonal elements of AR coefficient matrices.
    # If alpha!=0 and beta=NULL, beta will be estimated or alternately updated
    beta.oneshot = T, # whether estimate Gamma one-shot or alternately update Gamma
    converge = c("Xhat", "col.space"),
    thresh = 1e-05, # convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates
    maxit = 100, # maximum number of iterations
    ...
) {
  step2.b.update <- match.arg(step2.b.update)

  if (step1 && (!is.logical(use.ecm) || length(use.ecm) != 1 || is.na(use.ecm))) {
    stop("`use.ecm` must be a single TRUE/FALSE value in `siap`.")
  }
  if (!is.null(incl.cov)) {
    if (!is.logical(incl.cov) || length(incl.cov) != 1 || is.na(incl.cov)) {
      stop("`incl.cov` must be NULL or a single TRUE/FALSE value in `siap`.")
    }
    if (incl.cov != use.ecm) {
      warning("[siap] strict legacy mode: forcing `incl.cov <- use.ecm`.")
    }
  }
  incl.cov <- use.ecm

  # ----------- step 1 -----------
  if (step1 == T) {
    fit1 <- sia.prd_phase1(X, rank.max = rank.max, lambda = lambda, lambda_mu = lambda_mu1,
                           Z = step1.Z,
                           prd = step1.prd, use.ecm = use.ecm, incl.cov = incl.cov, ... )
    print("step 1 done.")
  } 
  
  # ----------- step 2 -----------
  if (step2 == T) {
    if (!is.null(fit1)) {
      x_imp1 = fit1$x_imp1
      if (is.null(x_imp1)) {
        x_imp1 <- fit1$a %*% t(fit1$b)
        if (isTRUE(fit1$prd) && !is.null(fit1$Theta) && !is.null(fit1$Phi)) {
          x_imp1 <- x_imp1 + t(fit1$Theta) %*% t(fit1$Phi)
        }
        x_imp1[!is.na(X)] <- X[!is.na(X)]
      }
      x_imp1 <- as.matrix(x_imp1)
      x_imp1[, apply(X, 2, function(v) all(is.na(v)))] <- NA
    } else {
      x_imp1 = X
    }
    
    if (step2.ar == F) {
      alpha = 0
    }
    fit2 <- sia.prd_phase2(x_imp1, rank.max = step2.rank, lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, lambda_mu = lambda_mu2, time.lag = time.lag,
                           Z = step2.Z, prd = step2.prd, b.update = step2.b.update, ...)
    print("step 2 done.")
  } else {
    fit2 <- NULL
  }
  
  return(list(fit1 = fit1, fit2 = fit2,
              step1 = step1, step2 = step2,
              step1.prd = step1.prd, use.ecm = use.ecm,
              step2.ar = step2.ar, step2.prd = step2.prd,
              step2.b.update = step2.b.update))
}

build_ssi_design <- function(
    dates,
    prd.type = c("l", "s", "y"),
    int.knots.ls = list(l = c(778.50, 1605.25, 2445.50)),
    phase_l = NULL,
    length.out = NULL
) {
  if (is.null(length.out)) {
    if (is.null(dates)) {
      stop("`dates` is required when `length.out` is NULL.")
    }
    length.out <- length(dates)
  }
  if (is.null(phase_l)) {
    if (is.null(dates) || length(dates) == 0 || is.na(dates[1])) {
      stop("`phase_l` is NULL and cannot be inferred because `dates[1]` is missing.")
    }
    phase_l <- lubridate::ymd(20191201) - lubridate::ymd(dates[1])
  }
  tx <- seq_len(length.out)
  get_Phi(tx = tx,
          type = prd.type,
          spline.type = "b",
          phase_l = phase_l,
          int.knots = list(l = int.knots.ls$l),
          length.out = length.out)
}

resolve_cov_controller <- function(incl.cov, use.ecm, context = "sia.prd_phase1") {
  if (!is.logical(incl.cov) || length(incl.cov) != 1 || is.na(incl.cov)) {
    stop("`incl.cov` must be a single TRUE/FALSE value.")
  }
  if (is.null(use.ecm)) {
    use.ecm <- incl.cov
  } else if (!is.logical(use.ecm) || length(use.ecm) != 1 || is.na(use.ecm)) {
    stop("`use.ecm` must be NULL or a single TRUE/FALSE value.")
  }
  if (!incl.cov && use.ecm) {
    warning(sprintf("[%s] `use.ecm=TRUE` requires covariance modeling; forcing `incl.cov=TRUE`.", context))
    incl.cov <- TRUE
  }
  list(incl.cov = incl.cov, use.ecm = use.ecm)
}

init_theta <- function(y, Phi, lambda_mu = 0) {
  y <- as.matrix(y)
  Phi <- as.matrix(Phi)
  if (nrow(y) != nrow(Phi)) {
    stop("`y` and `Phi` must have the same number of rows when initializing Theta0.")
  }
  q <- ncol(Phi)
  r <- ncol(y)
  Theta0 <- matrix(0, nrow = q, ncol = r)
  estimable <- rep(FALSE, r)

  for (j in seq_len(r)) {
    obs_j <- !is.na(y[, j])
    if (!any(obs_j)) {
      next
    }
    estimable[j] <- TRUE
    y_obs <- y[obs_j, j, drop = FALSE]
    Phi_obs <- Phi[obs_j, , drop = FALSE]
    if (lambda_mu == 0) {
      Theta0[, j] <- as.numeric(lm.fit(x = Phi_obs, y = y_obs)$coefficients)
    } else {
      gram <- crossprod(Phi_obs) + diag(lambda_mu, q)
      rhs <- crossprod(Phi_obs, y_obs)
      Theta0[, j] <- as.numeric(solve(gram, rhs))
    }
  }

  if (!any(estimable)) {
    stop("Cannot initialize Theta0: all response columns are fully missing.")
  }

  colnames(Theta0) <- colnames(y)
  rownames(Theta0) <- colnames(Phi)
  Theta0
}

sia.prd_phase2 <- function(
    X, # X only has chunk missingness
    rank.max,
    # --------- hyper-parameters ----------
    lambda1 = 0,
    lambda2 = 0,
    alpha = 0,
    lambda3 = 0,
    lambda_mu = 0,
    # ---------------- solving strategy ---------------
    beta = NULL, # AR parameters. A p x r matrix, indicating the diagonal elements of AR coefficient matrices.
    time.lag = 2, 
    beta.oneshot = T, # whether estimate Gamma or alternately update Gamma
    b.update = c("iterative", "whole", "auto"), # strategy for updating b when alpha != 0
    # If alpha!=0, beta.oneshot = F and beta=NULL, beta will be alternately updated
    Z = NULL, # optional design matrix for column-wise mean structure, d2 x q
    ## --------- splines ----------
    prd = T, # whether E(b_t) is a prd splines or constant 0
    prd.type = c("l", "s", "y"),
    int.knots.ls = list(l = c(778.50, 1605.25, 2445.50)),
    phase_l = NULL,
    # ----- A covariance (mostly discarded) ------
    incl.a_cov = F, # whether to inclue A covariance structure estimation
    a_cov = NULL, 
    glasso.rho = .2, # Regularization parameter for glasso.
                  # choose a rho in .01~.1 such that rho is small but large enough to give valid estimate.
    glasso.zero = NULL, # (Optional) indices of entries of inverse covariance to be constrained to be zero.
    glasso.approx = F, # Approximation flag: if true, computes Meinhausen-Buhlmann(2006) approximation
    # ------------- solving strategy -------------
    converge = c("Xhat", "col.space"),
    thresh = 1e-05, # convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates
    maxit = 100, # maximum number of iterations
    # -------- debug ---------
    fix.a = F,
    fix.b = F,
    fix.Theta = F,
    colrsdl.refit = F, # use provided column residuals to refit
    d_sigma = NULL,
    ...
  ) {
  
  converge = match.arg(converge)
  alg = "vanilla.prob"
  proj = F
  
  d1 = nrow(X)
  d2 = ncol(X)
  b.update <- match.arg(b.update)
  if (!is.numeric(lambda_mu) || length(lambda_mu) != 1 || is.na(lambda_mu) || lambda_mu < 0) {
    stop("`lambda_mu` must be a single non-negative numeric value.")
  }
  
  mean_mode <- "none"
  if (!is.null(Z)) {
    Phi <- as.matrix(Z)
    if (nrow(Phi) != d2) {
      stop("`Z` must have nrow(Z) == ncol(X).")
    }
    mean_mode <- "design"
  } else if (prd == T) {
    Phi <- build_ssi_design(
      dates = colnames(X),
      prd.type = prd.type,
      int.knots.ls = int.knots.ls,
      phase_l = phase_l,
      length.out = d2
    )
    mean_mode <- "legacy_prd"
  } else {
    Phi <- NULL
  }
  use_mean <- !is.null(Phi)
  
  if (is.null(beta)) {
    p = time.lag
  } else {
    p = nrow(beta) # beta: p x r
  }
  
  f.norm <- function(A, squared = F, na.rm = F){ # Frobenius norm of A
    if (squared == T) return(sum(A*A, na.rm = na.rm))
    return(sqrt(sum(A*A, na.rm = na.rm)))
  }

  solve_theta_two_sided <- function(R, A, Phi, lambda_mu = 0) {
    G_phi <- crossprod(Phi)
    G_a <- crossprod(A)
    svd_phi <- svd(G_phi)
    svd_a <- svd(G_a)
    mat_v = outer(svd_a$d, svd_phi$d) + lambda_mu # r x q
    RHS <- crossprod(Phi, t(R) %*% A) # q x r
    if (lambda_mu == 0) {
      return(chol2inv(chol(G_phi)) %*% RHS %*% chol2inv(chol(G_a)))
    }
    theta <- svd_phi$u %*% ((t(svd_phi$u) %*% RHS %*% svd_a$u) / t(mat_v)) %*% t(svd_a$u) # q x r
    return(theta)
  }
  
  Pi <- function(A, Omega = !is.na(X)) {
    # - Omega: default is the observed index set
    A[is.na(X)] = 0
    return(A*Omega)
  }
  
  pMat <- function(r, n){
    # permutation matrix generator
    # r can be any integer
    r = r %% n
    if (r==0) {
      pp = 1:n
    } else {
      pp = c((r+1):n, 1:r)
    }
    as(as.integer(pp), "pMatrix")
  }
  
  C0_func <- function(beta) {
    # C0(beta)
    C0 = Reduce("+", lapply(1:p, function(j) pMat(d2 - j, d2) %x% diag(beta[j,],ncol(beta)))) - pMat(d2, d2) %x% diag(rank.max)
    C0 = C0[(p*rank.max+1):(d2*rank.max), ]
    if (proj == T) {
      C0 = C0 * rep(Ome_p_1, each = rank.max)
    }
    return(C0)
  }
  
  M_func <- function(beta) {
    M = sqrt(alpha) * C0_func(beta)
  }
  
  Ja_func <- function() {
    if (alg == 'bayesian') {
      t = set2ind(union(comp(S2), 1:p))
      Ja = solve(pMat(d2, d2) %x% (t(a0)%*%a0) + lambda2 * (diag(t, d2) %x% diag(rank.max)))
    } else if (alg == 'separate') {
      t1 = set2ind(setdiff(union(comp(S2), 1:p), S1)) 
      t2 = set2ind(intersect(union(comp(S2), 1:p), S1)) 
      Ja = solve(pMat(d2, d2) %x% (t(a0)%*%a0) + 
                   lambda2 * (diag(t1, d2) %x% diag(rank.max)) +
                   lambda3 * (diag(t2, d2) %x% diag(rank.max)))
    } else if (alg == "vanilla.prob") {
      t = set2ind(1:p)
      # Ja = solve(pMat(d2, d2) %x% (t(a0)%*%a0) + lambda2 * (diag(t, d2) %x% diag(rank.max)))
      m_ls <- lapply(1:d2, function(jj) solve(t(a0)%*%a0 + as.numeric(jj<=p) * lambda2 * diag(rank.max)))
      Ja = Matrix::bdiag(m_ls)
    } else {
      Ja = pMat(d2, d2) %x% solve(t(a0)%*%a0 + lambda*diag(rank.max))
    }
    return(Ja)
  }
  
  get_gb <- function(t) {
    beta0 <- rbind(-1, beta)
    lam_k <- function(k,j) {
      if(k==j) {
        rep(0, ncol(b))
      } else {
        beta0[k+1,] * b[t+j-k,]
      }
    }
    contrast_j <- function(j) {
      rowSums(sapply(0:p, lam_k, j=j))
    }
    rowSums(sapply(max(0,p-t+1):min(p,d2-t), function(j) beta0[j+1,] * contrast_j(j)))
  }
  
  get_gg <- function(t) {
    beta0 <- rbind(-1, beta)
    rowSums(sapply(max(0,p-t+1):min(p,d2-t), function(j) beta0[j+1,]^2))
  }
  
  update_bt <- function(t) {
    solve(t(a0)%*%a0 + (t<=p)*diag(lambda2, nrow = ncol(b)) + alpha * diag(get_gg(t), nrow = ncol(b)))%*%(t(a0) %*% Xs[,t] - alpha * get_gb(t))
  }
  
  start_time = Sys.time()
  
  # ----------- initialization ----------------
  ## a0, b0
  si = softImpute::softImpute(X, rank.max = rank.max, lambda = lambda1)
  if (sum(si$d>0) < rank.max) {
    rank.max = min(rank.max, sum(si$d>0)) # corrected rank.max
    cat("Maximum rank changed to", rank.max, "by softImp.")
  }
  if (rank.max == 0) {
    stop("The chosen lambda is too large. Please decrease it.")
  }
  d0 = sqrt(si$d)[1:rank.max]
  u0 = as.matrix(as.matrix(si$u)[,1:rank.max])
  v0 = as.matrix(as.matrix(si$v)[,1:rank.max])
  a0 = t(t(u0) * d0)
  b0 = t(t(v0) * d0)
  
  get.prd <- function(type = c("last", "now"), transpose = F) {
    # return t(Theta0) %*% t(Phi) when a mean design is enabled
    type = match.arg(type)
    if (use_mean) {
      if (type == "last") {
        if (transpose == T) {
          return(Phi %*% Theta0)
        }
        return(t(Theta0) %*% t(Phi))
      } else {
        if (transpose == T) {
          return(Phi %*% Theta)
        }
        return(t(Theta) %*% t(Phi))
        }
    } else {
      if (transpose == T) {
        return(matrix(0, d2, rank.max))
      }
      return(matrix(0, rank.max, d2))
    }
  }
  
  if (use_mean) {
    ## Theta0
    S = which(apply(X, 2, function(x) all(is.na(x))))
    y <- b0
    y[S, ] <- NA # Do not use 0 value b's to get a better initialization
    Theta0 <- init_theta(y = y, Phi = Phi, lambda_mu = lambda_mu)
    Theta <- Theta0

    ## Re-Initialize b0
    #b0 = Phi  %*% Theta0
    b0[S,] = (Phi  %*% Theta0)[S,]
  } else {
    Theta0 = 0
    Theta = 0
  }
  
  
  # set functions ------------------
  comp <- function(S) {
    return(setdiff(1:d2, S))
  }
  
  set2ind <- function(S) {
    t = rep(0, d2)
    t[S] = 1
    return(t)
  }
  
  ind2set <- function(ind) {
    return(which(ind==1))
  }
  
  Ome_prime = apply(X, 2, function(x) all(is.na(x))) # The indicator set of columns that are entirely missing
  
  # estimate gamma -------------------
  if (use_mean) {
    b0_tilde = b0 - get.prd(transpose = T)
  } else {
    b0_tilde = b0
  }
  if (alpha != 0 && is.null(beta)) { # initialization of beta
    Ome_p_1 = Ome_prime # downtime left p neighbor indicator set. Length: n-p
    for (i in which(Ome_prime==1)){
      l = max(1, i-p)
      Ome_p_1[l:i] = 1
    }
    Ome_p_1 = Ome_p_1[1:(d2-p)]
    S1_trunc = which(Ome_p_1 == 1) # downtime left p neighbor index set, \subset [0,n-p]
    S1_comp_trunc = which(Ome_p_1 != 1) # S1 complement, \subset [0,n-p]
    if (length(S1_comp_trunc) == 0) {
      warning("No data points outside the lag window. Using non-downtime points within the lag window to estimate AR parameters.
      Consider decreasing `time.lag`.")
      Ome_p = Ome_prime[1:(d2-p)]
      S_comp_trunc = which(Ome_p == 0) # S: non-downtime columns, \subset [0, n-p]
      l = sapply(1:p, function(k) sapply(1:p, function(i) Reduce("+", lapply(S_comp_trunc, function(j) b0_tilde[j+p-i,] * b0_tilde[j+p-k,]))))
      L = array(l, dim = c(rank.max, p, p))
      R = sapply(1:p, function(k) Reduce("+", lapply(S_comp_trunc, function(j) as.matrix(b0_tilde[j+p,] * b0_tilde[j+p-k,]))))
    } else {
      l = sapply(1:p, function(k) sapply(1:p, function(i) Reduce("+", lapply(S1_comp_trunc, function(j) b0_tilde[j+p-i,] * b0_tilde[j+p-k,]))))
      L = array(l, dim = c(rank.max, p, p))
      R = sapply(1:p, function(k) Reduce("+", lapply(S1_comp_trunc, function(j) as.matrix(b0_tilde[j+p,] * b0_tilde[j+p-k,]))))
    }
    if (!is.array(R)) {
      R = matrix(R, nrow = 1)
    }
    beta = sapply(1:rank.max, function(i) solve(matrix(L[i,,], nrow = dim(L)[2], ncol = dim(L)[3]), R[i,]))
    beta = matrix(beta, nrow = time.lag, ncol = rank.max)
  }
  (rm(b0_tilde))
  
  # estimate A covariance matrix --------------------
  if (incl.a_cov == T) {
    # Use a0 and glasso to estimate a_cov
    s = a0 %*% t(a0) / ncol(a0)
    s_diag = diag(s)
    s_std = t(s/sqrt(s_diag))/sqrt(s_diag)
    cov.fit = glasso(s_std, rho = glasso.rho,
                     zero = glasso.zero, approx = glasso.approx) # choose a rho in .01~.1 such that rho is small but large enough to give valid estimate
    #a_cov = cov.fit$w
    a_cov_inv = t(cov.fit$wi*sqrt(s_diag))*sqrt(s_diag)
  }
  
  # loss function ----------
  lossFunc <- function() {
    
  } 
  
  for (i in 1:maxit) {
    X0 = a0 %*% t(b0)
    # Update Theta --------------------
    if (use_mean && fix.Theta == F) {
      Xs = Pi(X - a0 %*% t(b0)) + a0 %*% get.prd()
      # C0Phi = Reduce("+", lapply(1:p, function(j) (pMat(d2 - j, d2) %*% Phi) %x% diag(beta[j,],ncol(beta))))
      #           - (pMat(d2, d2) %*% Phi) %x% diag(rank.max)
      # C0Phi = C0Phi[(p*rank.max+1):(d2*rank.max), ] # C0 %*% Phi
      # 
      # mat1 = lambda2 * (t(Phi) %*% diag(set2ind(1:p)) %*% Phi) %x% diag(rank.max) + alpha * t(C0Phi) %*% C0Phi
      # mat2 = lambda2 * (t(Phi) %*% diag(set2ind(1:p))) %x% diag(rank.max) + 
      #           alpha * (t(C0Phi) %*% C0_func(beta))
      # Theta = chol2inv(chol(mat1)) %*% mat2 %*% as.vector(t(b0))
      # Theta = matrix(Theta, ncol = rank.max, nrow = 7, byrow = T)
      Theta = solve_theta_two_sided(R = Xs, A = a0, Phi = Phi, lambda_mu = lambda_mu)
    } else {
      Theta = Theta0
    }
    if (colrsdl.refit && use_mean) {
      # d_sigma
      Theta = solve_theta_two_sided(R = Xs, A = a0, Phi = Phi / d_sigma^2, lambda_mu = lambda_mu)
    }
    
    # Update b -----------------------
    if (fix.b == F) {
      Xs = Pi(X) - Pi(a0 %*% t(b0)) + a0 %*% t(b0) - a0 %*% get.prd("now")
      if (alpha != 0) {
        if (colrsdl.refit) {
          # d_sigma
          M = M_func(beta)
          m_ls <- lapply(1:d2, function(jj) solve(t(a0)%*%a0/d_sigma[jj]^2 + as.numeric(jj<=p) * lambda2 * diag(rank.max)))
          Ja = Matrix::bdiag(m_ls)
          b = (Ja - Ja %*% t(M) %*% solve(diag((d2-p)*rank.max) + M%*%Ja%*%t(M)) %*% M %*% Ja) %*% as.vector(t(a0) %*% t(t(Xs)/d_sigma^2))
        } else {
          update_b_mode <- b.update
          if (update_b_mode == "auto") {
            update_b_mode <- ifelse(rank.max <= 15, "whole", "iterative")
          }
          if (update_b_mode == "whole") { # update b as a whole
            C0 = C0_func(beta)
            M = M_func(beta)
            Ja = Ja_func() # J_A matrix
            
            b = (Ja - Ja %*% t(M) %*% solve(diag((d2-p)*rank.max) + M%*%Ja%*%t(M)) %*% M %*% Ja) %*% as.vector(t(a0) %*% Xs)
            
            b = matrix(b, ncol = rank.max, nrow = d2, byrow = T)
          } else { # update b1, b2, ... iteratively
            b = b0 - get.prd("now", transpose = T)
            for (j in 1:d2) {
              b[j,] <- update_bt(j) # b will be changed in the way
            }
          }
        }
      } else { # softImpute
        b = t(Xs) %*% a0 %*% solve(t(a0) %*% a0 + lambda2 * diag(1, rank.max))
      }
      
      b = b + get.prd("now", transpose = T)
      
    } else {
      b = b0
    }
    
    
    # Update a -----------------------
    if (fix.a == F) {
      Xs = t(Pi(X)) - t(Pi(a0 %*% t(b))) + b %*% t(a0)
      if (incl.a_cov == F) { # no covariance structure
        a = t(Xs) %*% b %*% solve(t(b) %*% b + lambda1 * diag(1, rank.max))
      } else { # Sigma_A
        a = chol2inv(chol((t(b) %*% b) %x% diag(d1) + lambda1 * diag(rank.max) %x% a_cov_inv)) %*% as.vector(t(Xs) %*% b)
        a = matrix(a, ncol = rank.max, nrow = d1)
      }
    } else {
      a = a0
    }
    if (colrsdl.refit) {
      a = t(Xs) %*% (b/d_sigma^2) %*% solve(t(b) %*% (b/d_sigma^2) + lambda1 * diag(1, rank.max))
    }
    
    # Update beta ------------------
    if (beta.oneshot == F) {
      if (use_mean) {
        b_tilde = b - get.prd(transpose = T)
      } else {
        b_tilde = b
      }
      l = sapply(1:p, function(k) sapply(1:p, function(i) Reduce("+", lapply(1:(d2-p), function(j) b_tilde[j+p-i,] * b_tilde[j+p-k,])) + lambda3 / alpha * (i == k)))
      L = array(l, dim = c(rank.max, p, p))
      R = sapply(1:p, function(k) Reduce("+", lapply(1:(d2-p), function(j) as.matrix(b_tilde[j+p,] * b_tilde[j+p-k,]))))
      if (!is.array(R)) {
        R = matrix(R, nrow = 1)
      }
      beta = sapply(1:rank.max, function(i) solve(L[i,,], R[i,]))
      beta = matrix(beta, nrow = time.lag, ncol = rank.max)
    }
    
    # convergence ------------------
    check = "Reached maximum iteration"
    if (converge == 'col.space' && sum((a%*%solve(t(a)%*%a)%*%t(a) - a0%*%solve(t(a0)%*%a0)%*%t(a0))^2) +
        sum((b%*%solve(t(b)%*%b)%*%t(b) - b0%*%solve(t(b0)%*%b0)%*%t(b0))^2) <= thresh) {
      check = "No change of column spaces"
      break
    }
    if (converge == 'Xhat' && f.norm(a%*%t(b) - X0, squared = T) / f.norm(X0, squared = T) <= thresh) {
      check = "No change of X hat"
      break
    }
    
    a0 = a
    b0 = b
    Theta0 = Theta
  }
  # ---------- end -----------
  runtime = Sys.time() - start_time
  
  x_imp = a %*% t(b)
  x_imp[which(!is.na(X))] = X[which(!is.na(X))]
  return(list(a=a, b=b, Theta = Theta, beta = beta, x_imp = x_imp, rank.max = rank.max, time.lag = time.lag, prd.type=prd.type,
              mean_mode = mean_mode, mean_design = Phi, mean_coef = if (use_mean) Theta else NULL,
              train.err = f.norm(X - a%*%t(b), na.rm = T)/sqrt(d1*d2), 
              iter=i, time=runtime, check=check, prd = prd,
              hyperparam = list(lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, lambda_mu = lambda_mu, b.update = b.update)))
}

sia.prd_phase1 <- function(
    X, # If there are NA's in X, X is considered partially observed
    rank.max,
    # --------- hyper-parameters ----------
    lambda = 0,
    lambda_mu = 0,
    # --------- variance-covariance ----------
    incl.cov = F, # If T, the covariance matrix is estimated/input and very likely not diagonal.
    use.ecm = NULL, # NULL means auto: use.ecm <- incl.cov
    cov.method = c("glasso", "diag+lowrank"),
    ## ------- diag+lowrank --------
    cov.multi.update = F,
    rankSig = NULL,
    ## ------- glasso -------
    cov.fit = NULL, # huge output
    glasso.rho = .2, # Regularization parameter for glasso.
    # choose a rho in .01~.1 such that rho is small but large enough to give valid estimate.
    glasso.zero = NULL, # (Optional) indices of entries of inverse covariance to be constrained to be zero.
    glasso.approx = F, # Approximation flag: if true, computes Meinhausen-Buhlmann(2006) approximation
    ## ---------------- splines ---------------
    Z = NULL, # optional design matrix for column-wise mean structure, d2 x q
    prd = T,
    prd.type = c("l", "s", "y"),
    int.knots.ls = list(l = c(778.50, 1605.25, 2445.50)),
    phase_l = NULL,
    # ---------------- termination ----------
    converge = c("Xhat", "col.space"),
    thresh = 1e-05, # convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates
    thresh.sig = 1e-3,
    maxit = 100, # maximum number of iterations
    # jitter = 1e-7,
    # --------- debug -------
    fix.a = F, fix.b = F, fix.Theta = F,
    ...
  ) {
  
  converge = match.arg(converge)
  cov.method = match.arg(cov.method)
  d1 = nrow(X)
  d2 = ncol(X)
  if (!is.numeric(lambda_mu) || length(lambda_mu) != 1 || is.na(lambda_mu) || lambda_mu < 0) {
    stop("`lambda_mu` must be a single non-negative numeric value.")
  }
  cov_ctl <- resolve_cov_controller(incl.cov = incl.cov, use.ecm = use.ecm, context = "sia.prd_phase1")
  incl.cov <- cov_ctl$incl.cov
  use.ecm <- cov_ctl$use.ecm
  
  mean_mode <- "none"
  if (!is.null(Z)) {
    Phi <- as.matrix(Z)
    if (nrow(Phi) != d2) {
      stop("`Z` must have nrow(Z) == ncol(X).")
    }
    mean_mode <- "design"
  } else if (prd == T) {
    Phi <- build_ssi_design(
      dates = colnames(X),
      prd.type = prd.type,
      int.knots.ls = int.knots.ls,
      phase_l = phase_l,
      length.out = d2
    )
    mean_mode <- "legacy_prd"
  } else {
    Phi <- NULL
  }
  use_mean <- !is.null(Phi)
  
  f.norm <- function(A, squared = F, na.rm = F){ # Frobenius norm of A
    if (squared == T) return(sum(A*A, na.rm = na.rm))
    return(sqrt(sum(A*A, na.rm = na.rm)))
  }
  
  Pi <- function(A, Omega = !is.na(X)) {
    # - Omega: default is the observed index set
    A[is.na(X)] = 0
    return(A*Omega)
  }
  
  
  # Initialization --------------------------
  start_time = Sys.time()
  
  if (use_mean) {
    ## ------------ Theta0
    ### Theta0
    y <- t(X)
    Theta0 <- init_theta(y = y, Phi = Phi, lambda_mu = lambda_mu)

    Theta = Theta0
  } else {
    # in order for the returned list to work
    Theta0 = 0
    Theta = 0
  } 
  
  get.prd <- function(type = c("last", "now")) {
    # return t(Theta) %*% t(Phi) when a mean design is enabled
    type = match.arg(type)
    if (use_mean) {
      return(t(Theta) %*% t(Phi))
    } else {
      return(0)
    }
  }
  
  ## -------- A0, B0
  
  si = softImpute::softImpute(X - get.prd(), rank.max = rank.max, lambda = lambda)
  
  if (sum(si$d>0) < rank.max) {
    rank.max = min(rank.max, sum(si$d>0)) # corrected rank.max
    cat("Maximum rank changed to", rank.max, "by softImp.")
  }
  if (rank.max == 0) {
    stop("The chosen lambda is too large. Please decrease it.")
  }
  d0 = sqrt(si$d)[1:rank.max]
  u0 = as.matrix(as.matrix(si$u)[,1:rank.max])
  v0 = as.matrix(as.matrix(si$v)[,1:rank.max])
  a0 = t(t(u0) * d0)
  b0 = t(t(v0) * d0)
  a = a0
  b = b0
  d = d0
  u = u0
  v = v0

  #X[which(is.na(X))] = 0 # fills NA's in X with 0's
  
  lossFunc <- function(X, a, b, beta, index = c("overall", "b", "a", "Theta")) { # Loss function
    index = match.arg(index)
    if (index == "overall") {
      return(f.norm)
    } else if (index == "b") {
      
    } else if (index == "a") {
      
    } else {
      
    }
  }
  
  # ----------- warm start ---------------
  for (i in 1:maxit) {
    Xhat = a %*% t(b)   # a == a0, b == b0 at iteration start
    prd  = get.prd()
    X0   = Xhat + prd

    if (fix.b == F) {
      Xs = Pi(X) - Pi(Xhat + prd) + Xhat
      b = t(Xs) %*% u %*% diag(d / (d ^ 2 + lambda), length(d))
      # compatibility
      svd.1 = svd(t(t(b)*d), nu = min(d1,rank.max), nv = min(d2,rank.max))
      v = svd.1$u
      d = sqrt(svd.1$d)
      b = t(t(v)*d)
      u = u %*% svd.1$v # compatibility a0%*%t(b)
      a = t(t(u)*d)
      Xhat = a %*% t(b)
    } else { # softImpute
      b = b0
      v = v0
      d = d0
    }

    if (fix.a == F) {
      Xs = t(Pi(X)) - t(Pi(Xhat + prd)) + t(Xhat)
      a = t(Xs) %*% v %*% diag(d / (d ^ 2 + lambda), length(d))
      # compatibility
      svd.1 = svd(t(t(a)*d), nu = min(d1,rank.max), nv = min(d2,rank.max))
      u = svd.1$u
      d = sqrt(svd.1$d)
      a = t(t(u)*d)
      v = v %*% svd.1$v # compatibility
      b = t(t(v)*d)
      Xhat = a %*% t(b)
    } else {
      a = a0
      u = u0
    }

    if (use_mean && fix.Theta == F) {
      Xs = Pi(X) - Pi(Xhat + prd) + prd
      Theta = solve(t(Phi) %*% Phi + diag(lambda_mu, ncol(Phi))) %*% t(Phi) %*% t(Xs)
    } else {
      Theta = Theta0
    }

    ## convergence -----------
    check = "Reached maximum iteration"
    if (converge == 'Xhat' && f.norm(Xhat + get.prd("now") - X0, squared = T) / f.norm(X0, squared = T) <= thresh) {
      check = "No change of X hat"
      break
    }
    
    v0 = v
    u0 = u
    d0 = d
    a0 = a
    b0 = b
    Theta0 = Theta
  }
  
  # -------- incl.cov==T, SigmaX -----------
  if (incl.cov==T) {
    ## ------- initialization ---------
    
    v0 = v
    u0 = u
    d0 = d
    a0 = a
    b0 = b
    Theta0 = Theta

    Omega = !is.na(X) # d1 x d2
    S_partial = which(apply(X, 2, function(x) any(!is.na(x)) & any(is.na(x)))) # partially missing
    S_full = which(apply(X, 2, function(x) all(!is.na(x)))) # fully observed
    S_mis_full <- which(apply(X,2,function(v) all(is.na(v)))) # fully missing
    
    if (cov.method == "glasso") {
      if (is.null(cov.fit)) {
        Sigma_inv = diag(1, d1)
        Sigma = diag(1, d1)
        # x.complete <- X[, S_full]
        # # nobs.x = Omega %*% t(Omega) # Number of pairwise-complete observations.
        # nobs.x <- ncol(x.complete)
        # mu.x = t(t(u0) * d0^2) %*% t(v0) + get.prd()
        # rsdl.x = X - mu.x
        # # rsdl.x[which(is.na(rsdl.x))] = 0
        # rsdl.x = rsdl.x[, S_full]
        # Sx = rsdl.x %*% t(rsdl.x) / nobs.x # pairwise-complete calculation can lead to non-positive covariance matrix
        # SDx = sqrt(diag(Sx))
        # Rx = t(Sx / SDx) / SDx
        # cor.fit = huge::huge(Rx, method = "glasso", cov.output = TRUE, lambda = glasso.rho, verbose = F) # where the most computation cost paid
        # CorX_inv = cor.fit$icov[1][[1]]
        # CorX = cor.fit$cov[1][[1]]
        # Sigma_inv = t(CorX_inv / SDx) / SDx
        # Sigma = t(cor.fit$cov[1][[1]] * SDx) * SDx
      } else {
        #CorX_inv = cor.fit$icov
        #CorX = cor.fit$cov
        Sigma_inv = cov.fit$icov 
        Sigma = cov.fit$cov 
      }
    } else {
      # diag + lowrank
      if (is.null(rankSig)) {
        rankSig = d1
      }
      Lambda.d <- c(rep(1, rankSig), rep(1, d1 - rankSig)) * 1e-4 
      # L <- rbind(diag(0.9, rankSig), matrix(0, d1 - rankSig, rankSig))
      set.seed(1)
      L <- matrix(rnorm(rankSig*d1), d1, rankSig) * 1e-2 
      Sigma <- diag(Lambda.d) + L%*%t(L)
      Sigma_inv <- diag(1/Lambda.d) - (L/Lambda.d)%*%solve(diag(1,rankSig)+t(L)%*%(L/Lambda.d))%*%t(L/Lambda.d)
    }
    
    

    # ------------ iteration --------------
    get_e_step <- function(j, get.b = F, second.order = T) {
      Omega_j = Omega[,j]
      Sigma_j = Sigma[Omega_j, Omega_j]
      if (F & sum(Omega_j)>d1/2 & sum(Omega_j)>1000) { # even slower. Discarded
        # if missing ratio of j-th column < 0.5
        Sigma_inv_arrange = Sigma_inv[c(which(Omega_j==F), which(Omega_j==T)), c(which(Omega_j==F), which(Omega_j==T))]
        m1 = sum(Omega_j==F) # number of missing entries
        Uj = matrix(0, nrow = d1, ncol = 2*m1)
        Uj[1:m1,1:m1] = diag(m1)
        Uj[(m1+1):d1,(m1+1):(2*m1)] = Sigma[which(Omega_j==T), which(Omega_j==F)]
        Vj = matrix(0, ncol = d1, nrow = 2*m1)
        Vj[(m1+1):(2*m1),1:m1] = diag(m1)
        Vj[1:m1,(m1+1):d1] = Sigma[which(Omega_j==F), which(Omega_j==T)]
        Sigma_j_inv = Sigma_inv_arrange + (Sigma_inv_arrange%*%Uj)%*% solve(diag(2*m1)-Vj%*%Sigma_inv_arrange%*%Uj) %*%(Vj%*%Sigma_inv_arrange)
        Sigma_j_inv = Sigma_j_inv[(m1+1):d1,(m1+1):d1]
      } else {
        if (cov.method == "diag+lowrank") {
          Lambda_j <- Lambda.d[Omega_j]
          L_j <- L[Omega_j,]
          Sigma_j_inv <- diag(1/Lambda_j, length(Lambda_j)) - (L_j/Lambda_j)%*%solve(diag(1,rankSig)+t(L_j)%*%(L_j/Lambda_j))%*%t(L_j/Lambda_j)
        } else {
          tryerr <- try({Sigma_j_inv = chol2inv(chol(Sigma_j))}, silent = TRUE) # 3s
          # if (inherits(tryerr, "try-error")) {
          #   Sigma_j_inv = chol2inv(chol(Sigma_j + diag(jitter, nrow(Sigma_j))))
          # }
        }
      }
      
      a_j = a[Omega_j, , drop = F]
      x_j_obs = X[Omega_j, j]
      if (use_mean) {
        Theta_j = Theta[, Omega_j, drop=F]
        if (get.b == T) {
          b_j = chol2inv(chol(t(a_j)%*%Sigma_j_inv%*% a_j + lambda*diag(rank.max))) %*%t(a_j)%*%Sigma_j_inv %*%(x_j_obs - t(Theta_j)%*%Phi[j,])
          b[j, ] <<- b_j
        }
        
        if (use.ecm == T) {
          # first order
          mu_obs = a[Omega_j, , drop = F] %*% b[j,] + t(Theta_j)%*%Phi[j,]
          mu_mis = a[!Omega_j, , drop = F] %*% b[j,] + t(Theta[, !Omega_j])%*%Phi[j,]
          x_j_mis <- mu_mis + Sigma[!Omega_j, Omega_j, drop=F] %*% Sigma_j_inv %*% (x_j_obs - mu_obs)
          #x_j <- X[, j]
          #x_j[!Omega_j, j] <- x_j_mis
          x_imp1[!Omega_j,j] <<- x_j_mis
        }
      } else {
        if (get.b == T) {
          b_j = chol2inv(chol(t(a_j)%*%Sigma_j_inv%*% a_j + lambda*diag(rank.max))) %*%t(a_j)%*%Sigma_j_inv %*%(x_j_obs)
          b[j, ] <<- b_j
        }
        
        if (use.ecm == T) {
          # first order
          mu_obs = a[Omega_j, , drop = F] %*% b[j,]
          mu_mis = a[!Omega_j, , drop = F] %*% b[j,]
          x_j_mis <- mu_mis + Sigma[!Omega_j, Omega_j, drop=F] %*% Sigma_j_inv %*% (x_j_obs - mu_obs)
          #x_j <- X[, j]
          #x_j[!Omega_j, j] <- x_j_mis
          x_imp1[!Omega_j,j] <<- x_j_mis
        } 
      }
      
      if (second.order == T) {
        V_j_mis <- Sigma[!Omega_j, !Omega_j, drop=F] - Sigma[!Omega_j, Omega_j, drop=F] %*% Sigma_j_inv %*% Sigma[Omega_j, !Omega_j, drop=F]
        xx_j_mis <- V_j_mis + x_j_mis %*% t(x_j_mis)
        xx_j <- x_imp1[,j] %*% t(x_imp1[,j])
        xx_j[!Omega_j, !Omega_j] <- xx_j_mis
        return(xx_j)
      } 
    }
    
    for (i in 1:maxit) {
      # psi = sqrt(diag(Sigma))
      X0 = a0 %*% t(b0) + get.prd()
      Xhat <- a %*% t(b) + get.prd()
      if (use.ecm == T) {
        # ------- E step -------
        x_imp1 <- X
        # XX_imp1 <- Reduce("+", lapply(1:d2, get_e_step, get.b = F))
        XX_imp1 <- Reduce(function(acc, j) acc + get_e_step(j, get.b = F), S_partial, init = matrix(0, d1, d1))
        if (is.null(XX_imp1)) XX_imp1 <- matrix(0, nrow = d1, ncol = d1)
        x_imp1[, S_mis_full] <- Xhat[, S_mis_full]
        XX_imp1 <- XX_imp1 + length(S_mis_full) * Sigma + (a %*% t(b) + get.prd())[, S_mis_full] %*% t((a %*% t(b) + get.prd())[, S_mis_full])
        XX_imp1 <- XX_imp1 + x_imp1[, S_full] %*% t(x_imp1[, S_full])
        
        Sx = XX_imp1 + tcrossprod(Xhat) - crossprod(t(Xhat), t(x_imp1)) - crossprod(t(x_imp1), t(Xhat))
      }
      
      if (fix.b == F) {
        if (use.ecm == T) {
            b = t(x_imp1 - get.prd()) %*% Sigma_inv %*% a %*% solve(t(a) %*% Sigma_inv %*% a + diag(lambda, rank.max))
        } else {
          if (length(S_full) > 0) {
            if (use_mean) {
              b[S_full,] <- t(chol2inv(chol(t(a) %*% Sigma_inv %*% a + lambda*diag(rank.max)))%*%
                                t(a) %*% Sigma_inv %*% (X[,S_full]-t(Theta)%*%t(Phi[S_full,,drop=F])))
            } else {
              b[S_full,] <- t(chol2inv(chol(t(a) %*% Sigma_inv %*% a + lambda*diag(rank.max)))%*%
                                t(a) %*% Sigma_inv %*% (X[,S_full]))
            }
          }
          if (length(S_partial) > 0) {
            XX_imp1 <- Reduce(function(acc, j) acc + get_e_step(j, get.b = T), S_partial, init = matrix(0, d1, d1))
          }
          x_imp1[, S_mis_full] <- (a %*% t(b) + get.prd())[, S_mis_full]
          XX_imp1 <- XX_imp1 + length(S_mis_full) * Sigma + (a %*% t(b) + get.prd())[, S_mis_full] %*% t((a %*% t(b) + get.prd())[, S_mis_full])
        }
        
        # compatibility
        svd.1 = svd(t(t(b)*d0), nu = min(d1,rank.max), nv = min(d2,rank.max))
        v = svd.1$u
        d = sqrt(svd.1$d)
        b = t(t(v)*d)
        u = u %*% svd.1$v # compatibility a%*%t(b)
        a = t(t(u0)*d)
      } else { # softImpute
        b = b0
        v = v0
        d = d0
      }
      
      ### ---- update Sigma -----
      if (incl.cov == T){
        if (cov.method == "diag+lowrank") {
          while (1) {
            # an inner mini-ECM iteration
            ## E step
            Sig_z <- solve(diag(1, rankSig) + t(L) %*% (L / Lambda.d))
            xz <- Sx %*% (L/Lambda.d) %*% Sig_z
            
            ## M step
            L0 <- L
            L <- xz %*% solve(d2 * Sig_z + Sig_z%*%t(L/Lambda.d)%*% Sx %*%(L/Lambda.d)%*%Sig_z)
            
            Lambda.d0 <- Lambda.d
            Lambda.d <- diag(Sx - L%*%t(xz)) / d2
            #print(paste("multi.update",i, f.norm(Lambda.d0 - Lambda.d, squared = T) / f.norm(Lambda.d, squared = T), f.norm(L0 - L, squared = T) / f.norm(L, squared = T)))
            
            if (cov.multi.update == F || (f.norm(Lambda.d0 - Lambda.d, squared = T) / f.norm(Lambda.d, squared = T) <= thresh.sig &&f.norm(L0 - L, squared = T) / f.norm(L, squared = T) <= thresh.sig)) break
          }
          
          
          Sigma0 <- diag(Lambda.d0, length(Lambda.d0)) + L0%*%t(L0) # last iteration
          Sigma <- diag(Lambda.d, length(Lambda.d)) + L%*%t(L)
          Sigma_inv <- diag(1/Lambda.d, length(Lambda.d)) - (L/Lambda.d)%*%solve(diag(1,rankSig)+t(L)%*%(L/Lambda.d))%*%t(L/Lambda.d)
          
        } else {
          # glasso
          Sigma0 <- Sigma
          Sx = XX_imp1 + (a %*% t(b) + get.prd()) %*% t(a %*% t(b) + get.prd()) -
            (a %*% t(b) + get.prd()) %*% t(x_imp1) - x_imp1 %*% t(a %*% t(b) + get.prd())  # ecm
          Sx <- Sx / d2
          SDx = sqrt(diag(Sx))
          Rx = t(Sx / SDx) / SDx
          cor.fit = huge::huge(Rx, method = "glasso", scr = T, cov.output = TRUE, lambda = glasso.rho, verbose = F) # where the most computation cost paid
          CorX_inv = cor.fit$icov[1][[1]]
          CorX = cor.fit$cov[1][[1]]
          Sigma_inv = t(CorX_inv / SDx) / SDx
          Sigma = t(cor.fit$cov[1][[1]] * SDx) * SDx
        }
      }
      
      if (fix.a == F) {
        # Xs = Pi(X / psi) - Pi((a %*% t(b) + get.prd()) / psi) + a %*% t(b) / psi
        # a = Xs %*% b / outer(lambda*psi^2, d^2, FUN = "+")
        # a = a * psi
        a = (x_imp1 - get.prd()) %*% b %*% solve(t(b) %*% b + diag(lambda, rank.max))
        
        # compatibility
        svd.1 = svd(t(t(a)*d), nu = min(d1,rank.max), nv = min(d2,rank.max))
        u = svd.1$u
        d = sqrt(svd.1$d)
        a = t(t(u)*d)
        v = v %*% svd.1$v # compatibility
        b = t(t(v)*d)
      } else {
        a = a0
        u = u0
      }
      
      if (use_mean && fix.Theta == F) {
        # Xs = Pi(X / psi) - Pi((a %*% t(b) + get.prd()) / psi) + get.prd() / psi
        # Theta = chol2inv(chol(t(Phi) %*% Phi)) %*% t(Phi) %*% t(Xs)
        # Theta = t(t(Theta) * psi)
        Theta = chol2inv(chol(t(Phi) %*% Phi + diag(lambda_mu, ncol(Phi)))) %*% t(Phi) %*% t(x_imp1 - a %*% t(b))
        
      } else {
        Theta = Theta0
      }
      
      
      
      #print(f.norm(X  - a%*%t(b) - t(Theta) %*% t(Phi), na.rm = T))
      check = "Reached maximum iteration"
      if (converge == 'col.space' && sum((a%*%solve(t(a)%*%a)%*%t(a) - a0%*%solve(t(a0)%*%a0)%*%t(a0))^2) +
          sum((b%*%solve(t(b)%*%b)%*%t(b) - b0%*%solve(t(b0)%*%b0)%*%t(b0))^2) <= thresh) {
        check = "No change of column spaces"
        break
      }
      if (converge == 'Xhat' && f.norm(a%*%t(b) + get.prd("now") - X0, squared = T) / f.norm(X0, squared = T) <= thresh) {
        check = "No change of X hat"
        
        if (cov.method == "diag+lowrank") {
          print(paste(i, f.norm(Sigma0 - Sigma, squared = T) / f.norm(Sigma0, squared = T),
                      f.norm(Lambda.d0 - Lambda.d, squared = T) / f.norm(Lambda.d, squared = T),
                      f.norm(L0 - L, squared = T) / f.norm(L, squared = T)))
          if (f.norm(Sigma0 - Sigma, squared = T) / f.norm(Sigma0, squared = T) <= thresh.sig &&
              f.norm(Lambda.d0 - Lambda.d, squared = T) / f.norm(Lambda.d, squared = T) <= thresh.sig) {
            check = "No change of X hat and Sigma hat"
            break
          }
        } else {
          print(paste(i, f.norm(Sigma0 - Sigma, squared = T) / f.norm(Sigma0, squared = T)))
          if (f.norm(Sigma0 - Sigma, squared = T) / f.norm(Sigma0, squared = T) <= thresh.sig) {
            check = "No change of X hat and Sigma hat"
            break
          }
        }
      }
      
      v0 = v
      u0 = u
      d0 = d
      a0 = a
      b0 = b
      Theta0 = Theta
    }
  }
  
  # ---------- final ecm ---------
  if (use.ecm == T) { 
    sapply(S_partial, get_e_step, get.b = F, second.order = F) # final update of x_imp1
  }
  
  # ------------ end -----------
  runtime = Sys.time() - start_time
  
  if (incl.cov == F) {
    return(list(a=a, b=b, u=u, v=v, d=d, Theta = Theta, Phi=Phi, rank.max = rank.max,
                mean_mode = mean_mode, mean_design = Phi, mean_coef = if (use_mean) Theta else NULL,
                train.err = f.norm(X  - a%*%t(b) - get.prd("now"), na.rm = T)/sqrt(d1*d2), 
                iter=i, time=runtime, check=check, prd=prd, lambda = lambda, lambda_mu = lambda_mu, prd.type = prd.type))
  } else if (incl.cov == T && use.ecm == T) {
    if (cov.method == "glasso") {
      return(list(a=a, b=b, u=u, v=v, d=d, Theta = Theta, Phi=Phi, rank.max = rank.max, Sigma = Sigma, x_imp1 = x_imp1, #x_mis = x_mis,
                  mean_mode = mean_mode, mean_design = Phi, mean_coef = if (use_mean) Theta else NULL,
                  train.err = f.norm(X  - a%*%t(b) - get.prd("now"), na.rm = T)/sqrt(d1*d2), glasso.rho = glasso.rho,
                  iter=i, time=runtime, check=check, prd=prd, lambda = lambda, lambda_mu = lambda_mu, prd.type = prd.type))
    } else {
      return(list(a=a, b=b, u=u, v=v, d=d, Theta = Theta, Phi=Phi, rank.max = rank.max, Sigma = Sigma, x_imp1 = x_imp1, #x_mis = x_mis,
                  mean_mode = mean_mode, mean_design = Phi, mean_coef = if (use_mean) Theta else NULL,
                  train.err = f.norm(X  - a%*%t(b) - get.prd("now"), na.rm = T)/sqrt(d1*d2), rankSig = rankSig,
                  iter=i, time=runtime, check=check, prd=prd, lambda = lambda, lambda_mu = lambda_mu, prd.type = prd.type))
    }
    
  } else {
    return(list(a=a, b=b, u=u, v=v, d=d, Theta = Theta, Phi=Phi, rank.max = rank.max, Sigma = Sigma, 
                mean_mode = mean_mode, mean_design = Phi, mean_coef = if (use_mean) Theta else NULL,
                train.err = f.norm(X  - a%*%t(b) - get.prd("now"), na.rm = T)/sqrt(d1*d2), 
                iter=i, time=runtime, check=check, prd=prd, lambda = lambda, lambda_mu = lambda_mu, prd.type = prd.type))
  }
  
}

get_Phi <- function(n.knots = c(l = 3, s = 2, y = 3),
                    phase_l = 0, # the min or max point in the period. Need to be numeric
                    type = c("l", "s", "y"), # control which periodicity to include
                    spline.type = c("m", "b"),
                    integrate = F,
                    length.out = 11 * 365,
                    tx = seq(1, length.out), # Need to be numeric
                    knots.selection = c("equidistant", "specify.l"),
                    int.knots = NULL # a list specifying knots
) {
  spline.type = match.arg(spline.type)
  knots.selection = match.arg(knots.selection)
  
  period = c(l=11*(4*365+1)/4, s=26.8, y=(4*365+1)/4)
  
  if (knots.selection == "equidistant") {
    int.knots_l = seq(0, period[1], length.out = n.knots["l"] + 2)[1+1:n.knots["l"]] 
    int.knots_s = seq(0, period[2], length.out = n.knots["s"] + 2)[1+1:n.knots["s"]]
    int.knots_y = seq(0, period[3], length.out = n.knots["y"] + 2)[1+1:n.knots["y"]]
  } else if (knots.selection == "specify.l"){ # specify long-term knots
    int.knots_l = int.knots$l
    int.knots_s = seq(0, period[2], length.out = n.knots["s"] + 2)[1+1:n.knots["s"]]
    int.knots_y = seq(0, period[3], length.out = n.knots["y"] + 2)[1+1:n.knots["y"]]
  }
  
  int.knots = list(
    l = int.knots_l,
    s = int.knots_s,
    y = int.knots_y
  )
  
  if (spline.type == "m") {
    # normalized splines
    intercept = T
    Phi = NULL
    for (i in type) {
      if (i == "l") {
        
        Phi = cbind(Phi, splines2::mSpline(x = tx, 
                                             periodic = TRUE, 
                                             degree = 3, 
                                             intercept = intercept,
                                             knots = as.numeric(phase_l)%%period[i] + int.knots[[i]],
                                             Boundary.knots = as.numeric(phase_l)%%period[i] + c(0, period[i])))
      } else {
        Phi = cbind(Phi, splines2::mSpline(x = tx, 
                                             periodic = TRUE, 
                                             degree = 3, 
                                             intercept = intercept,
                                             knots = int.knots[[i]],
                                             Boundary.knots = c(0, period[i])))
      }
        intercept = F # no intercept after the first iter
    }
    
  } else {
    # un-normalized splines
    intercept = T
    Phi = NULL
    for (i in type) {
      if (i == "l") {
        Phi = cbind(Phi, splines2::bSpline(x = tx, 
                                           periodic = TRUE, 
                                           degree = 3, 
                                           intercept = intercept,
                                           knots = as.numeric(phase_l)%%period[i] + int.knots[[i]],
                                           Boundary.knots = as.numeric(phase_l)%%period[i] + c(0, period[i])))
      } else {
        Phi = cbind(Phi, splines2::bSpline(x = tx, 
                                           periodic = TRUE, 
                                           degree = 3, 
                                           intercept = intercept,
                                           knots = int.knots[[i]],
                                           Boundary.knots = c(0, period[i])))
      }
      intercept = F # no intercept after the first iter
    }
  }
  if (integrate == T) {
    iPhi = update(Phi, integral = T)
    return(list(iPhi = iPhi, Phi = Phi))
  }
  return(Phi)
}
