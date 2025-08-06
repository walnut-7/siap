marss.wrapper = function(data, job, instance, model, p, r,...) {
  source("code/preprocess_functions.R")
  source("code/postprocess_functions.R")
  source("code/marss_model_specification.R")
  x <- data
  xt <- data[instance$div, , drop = F]
  x[instance$miss] <- NA
  x_big <- x
  x <- x[instance$div, , drop = F]
  pre_par <- preprocess_func(x, box.cox = F)
  x_work <- pre_par$x_ready
  d1 <- nrow(pre_par$x_ready)
  d2 <- ncol(pre_par$x_ready)

  ## -------- model specification ----------
  #model.list <- marss_model(d1 = d1, p = p, model = "simple")
  model.list <- marss_model(d1 = d1, p = p, r = r, model = model)
  
  ## -------------- fit model --------------
  start <- Sys.time()
  system.time({
    fit0 <- MARSS::MARSS(x_work, model=model.list, control=list(maxit=15), silent=T) # kem initialize. minit = 15.
    fit <- MARSS::MARSS(x_work, model=model.list, method="BFGS", inits = fit0, silent=T)
  })
  runtime.marss <- Sys.time() - start
  
  x_pred <- fit$ytT
  x_pred <- postprocess_func(x_pred, pre_par, x)
  
  start <- Sys.time()
  result <- try(boot.fit <- MARSS::MARSSboot(fit, nboot = 30), silent = T)
  runtime.boot <- Sys.time() - start
  
  if (inherits(result, "try-error")) {
    print("Bootstrap estimation error.")
    start <- Sys.time()
    result <- try(boot.fit <- MARSS::MARSSboot(fit, nboot = 20), silent = T)
    runtime.boot <- Sys.time() - start
  } 
  par_names <- matrix(unlist(strsplit(rownames(boot.fit$boot.params), split = '[.]')), ncol = 2, byrow = T)
  fit_new <- fit
  bootsamples <- array(dim = c(d1, d2, boot.fit$nboot))
  x_boots <- array(dim = c(nrow(x), ncol(x), boot.fit$nboot))
  sum_squares <- matrix(0, nrow(x), ncol(x))
  for (i in 1:(boot.fit$nboot)) {
    par_new <- fit$par
    for (j in 1:nrow(par_names)) {
      par_new[[par_names[j, 1]]][par_names[j, 2], ] <- boot.fit$boot.params[j, i]
    }
    fit_new$par <- par_new
    fit.hatyt_new <- MARSS::MARSShatyt(fit_new)
    bootsamples[,,i] <- fit.hatyt_new$ytT
    x_boot <- postprocess_func(bootsamples[,,i], pre_par, x)
    x_boots[,,i] <- x_boot
    sum_squares <- sum_squares + (x_pred - x_boot)^2
  }
  
  x_sd <- sqrt(1/(dim(bootsamples)[3])*sum_squares)
  
  ## ----------- return ------------  
  S.test.o <- which(is.na(x)&(!is.na(xt)), arr.ind = T) # S.test.o are all observed in data
  S.test.w <- setdiff(which(apply(x_big, 2, function(v) all(is.na(v)))), which(apply(data, 2, function(v) all(is.na(v))))) # may contain scattered missing entries in data
  S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
  
  ttl_coverage.w <- sum(abs(x_pred[, S.test.w] - xt[, S.test.w]) <= qnorm(0.975)*x_sd[,S.test.w], na.rm = T)
  ttl_coverage.o <- sum(abs(x_pred[S.test.o] - xt[S.test.o]) <= qnorm(0.975)*x_sd[S.test.o])
  
  ave_coverage_rate <- (ttl_coverage.w + ttl_coverage.o)/sum(is.na(x)&(!is.na(xt)))
  
  rr <- abs((x_pred - xt)/xt)
  rae.w <- rr[,S.test.w,drop=F]
  
  mask.o <- matrix(NA, nrow(x), ncol(x))
  mask.o[S.test.o] <- 1
  rae.o <- rr*mask.o
  
  list(x = x, x_pred=x_pred, x_sd=x_sd, x_boots = x_boots,
       rae = list(rae.w=rae.w, rae.o=rae.o), S.test.w = S.test.w, S.test.o = S.test.o,
       ave_coverage = ave_coverage_rate,
       ttl_coverage.w = ttl_coverage.w, n.w = (sum(is.na(x)&(!is.na(xt)))-nrow(S.test.o)),
       ttl_coverage.o = ttl_coverage.o, n.o = nrow(S.test.o),
       runtime = list(runtime.marss=runtime.marss, runtime.boot=runtime.boot),
       ytT=fit$ytT, bootsamples=bootsamples, instance = instance)
}
