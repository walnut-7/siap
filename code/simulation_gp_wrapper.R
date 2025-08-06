gp.wrapper <- function(data, job, instance, covfun, ...) { 
  library(GpGp)
  source("code/preprocess_functions.R")
  source("code/postprocess_functions.R")
  
  x <- data
  x[instance] <- NA
  pre_par <- preprocess_func(x, box.cox = F)
  d1 <- nrow(x)
  d2 <- ncol(x)
  
  x_df <- data.frame(
    row = rep(1:nrow(pre_par$x_ready), times = ncol(pre_par$x_ready)),
    col = rep(1:ncol(pre_par$x_ready), each = nrow(pre_par$x_ready)),
    value = as.vector(pre_par$x_ready)
  )
  
  not_missing <- !is.na(x_df$value)
  x_obs <- x_df$value[not_missing]
  locs_obs <- cbind(x_df$row[not_missing], x_df$col[not_missing])
  
  
  period = c(11 * 365, 26.8)
  tx = seq(1, ncol(pre_par$x_ready))
  Phi = cbind(
    splines2::mSpline(x = tx, df = 4, 
                      periodic = TRUE, 
                      degree = 3, 
                      intercept = T,
                      Boundary.knots = c(0, period[1])), # Phi_l
    splines2::mSpline(x = tx, df = 3, 
                      periodic = TRUE, 
                      degree = 3, 
                      intercept = F,
                      Boundary.knots = c(0, period[2])) # Phi_s
  ) # n x 7
  z <- cbind(Phi[locs_obs[, 2], ], locs_obs[, 1]) # n_obs x 8
  
  
  start <- Sys.time()
  gpfit <- GpGp::fit_model(y = x_obs, locs = locs_obs, X = z, covfun_name = covfun, m_seq = c(10,30))
  runtime <- Sys.time() - start
  
  is_missing <- is.na(x_df$value)
  locs_pred <- cbind(x_df$row[is_missing], x_df$col[is_missing])
  z_pred <- cbind(Phi[locs_pred[, 2], ], locs_pred[, 1]) # n_obs x 8
  
  start <- Sys.time()
  pred <- GpGp::predictions(fit = gpfit, locs_pred = locs_pred, X_pred = z_pred, m = 30) # m: number of neighbors
  runtime.pred <- Sys.time() - start
  
  obs_and_pred <- x_df$value
  obs_and_pred[is_missing] <- pred
  
  ncondsim <- 30
  start <- Sys.time()
  sims <- GpGp::cond_sim(fit = gpfit, locs_pred = locs_pred, X_pred = z_pred, nsims = ncondsim, m = 30)
  runtime.sims <- Sys.time() - start
  
  x_pred <- array(obs_and_pred, c(nrow(pre_par$x_ready), ncol(pre_par$x_ready)))
  x_pred <- postprocess_func(x_pred, pre_par, x)
  sum_squares <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  x_sims <- array(dim = c(nrow(x), ncol(x), ncondsim))
  for (j in 1:ncondsim) {
    obs_and_sim <-x_df$value
    obs_and_sim[is_missing] <- sims[, j]
    x_sim <- array(obs_and_sim, c(nrow(pre_par$x_ready), ncol(pre_par$x_ready)))
    x_sim <- postprocess_func(x_sim, pre_par, x)
    x_sims[,,j] <- x_sim
    sum_squares <- sum_squares + (x_pred - x_sim)^2
  }
  x_sd <- sqrt(1/ncondsim*sum_squares)
  
  
  ## ----------- return ------------  
  S.test.o <- which(is.na(x), arr.ind = T)
  S.test.w <- which(apply(x, 2, function(x) all(is.na(x))))
  S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
  
  ttl_coverage.w <- sum(abs(x_pred[, S.test.w] - data[, S.test.w]) <= qnorm(0.975)*x_sd[,S.test.w])
  ttl_coverage.o <- sum(abs(x_pred[S.test.o] - data[S.test.o]) <= qnorm(0.975)*x_sd[S.test.o])
  
  ave_coverage_rate <- (ttl_coverage.w + ttl_coverage.o)/sum(is.na(x))
  
  ### ----------- spectral mrae-------------
  rr <- (x_pred - data)/data
  mrae.w <- apply(rr[,S.test.w], 1, function(v) mean(abs(v)))
  
  mask.o <- matrix(NA, d1, d2)
  mask.o[S.test.o] <- 1
  mrae.o <- apply(rr*mask.o, 1, function(v) mean(abs(v), na.rm=T))
  
  list(x = x, x_pred = x_pred, #x_sims = x_sims, 
       x_sd = x_sd, mrae = list(mrae.w=mrae.w, mrae.o=mrae.o),
       ave_coverage = ave_coverage_rate, ave_coverage.w = ttl_coverage.w/(d1*length(S.test.w)), ave_coverage.o = ttl_coverage.o/nrow(S.test.o),
       runtime = list(runtime.gp=runtime, runtime.pred=runtime.pred, runtime.sims=runtime.sims),
       # gp = gpfit, pred = pred
       )
}
