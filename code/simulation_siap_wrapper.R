siap.wrapper <- function(data, job, instance, cp.alpha, preprocess = F,...) { 
  library(doParallel)
  library(Matrix)
  source("code/siap.R")
  source("code/cp.R")
  source("code/preprocess_functions.R")
  source("code/postprocess_functions.R")
  x <- data
  x[instance$miss] <- NA
  d1 <- nrow(x)
  d2 <- ncol(x)
  
  ## -------------- fit model --------------
  step1.prd <- list(...)$step1.prd
  step2.prd <- list(...)$step2.prd
#  if (is.null(preprocess)) {
#    preprocess <- ifelse(step1.prd==F | step2.prd==F, T, F)
#  }
  
  S.cal.w <- instance$cp$S.cal.w
  S.cal.o <- instance$cp$S.cal.o
  
  M.o <- !is.na(x)
  rownames(M.o) <- 1:nrow(x)
  M.o[, S.cal.w] <- F
  M.o[S.cal.o] <- F
  
  x.tr <- x
  x.tr[which(M.o == F)] <- NA

  if (preprocess) {
    pre_par <- preprocess_func(x.tr, diff = F, box.cox = F)
    x.tr <- pre_par$x_ready
  } else {
    pre_par <- NULL
  }
  
  out <- siap(x.tr, ...)
  
  x_imp = out$fit2$x_imp
  if (is.null(x_imp)) {
    x_imp = out$fit1$x_imp1
    if (is.null(x_imp)) {
      if (out$fit1$prd == T) {
        x_imp = out$fit1$a %*% t(out$fit1$b)+ t(out$fit1$Theta) %*% t(out$fit1$Phi)
      } else {
        x_imp = out$fit1$a %*% t(out$fit1$b)
      }
      x_imp[!is.na(x.tr)] = x.tr[!is.na(x.tr)]
    }
  }
  
  if (preprocess) {
    x_imp <- postprocess_func(x_imp, pre_par)
  }
  
  # ----------- cp ------------
  cp_res <- compute_cp(
    x = x,
    x_imp = x_imp,
    xt = data,
    cp.alpha = cp.alpha,
    S.cal.w = S.cal.w,
    S.cal.o = S.cal.o
  )
  S.test.o <- cp_res$S.test.o
  S.test.w <- cp_res$S.test.w
  ttl_coverage.w <- cp_res$ttl_coverage.w
  ttl_coverage.o <- cp_res$ttl_coverage.o
  ave_coverage_rate <- cp_res$ave_coverage
  cp_ls <- cp_res$cp
  
  # ----------- spectral mrae-------------
  #  calculated on test set but not S.test.o2w
  rr <- (x_imp - data)/data
  mrae.w <- apply(rr[,S.test.w], 1, function(v) mean(abs(v)))
  
  mask.o <- matrix(NA, d1, d2)
  mask.o[S.test.o] <- 1
  mrae.o <- apply(rr*mask.o, 1, function(v) mean(abs(v), na.rm=T))
  
  list(#x = x0, 
       mrae = list(mrae.w=mrae.w, mrae.o=mrae.o), 
       ave_coverage = ave_coverage_rate, 
       ave_coverage.w = ttl_coverage.w/(sum(is.na(x))-nrow(S.test.o)), 
       ave_coverage.o = ttl_coverage.o/nrow(S.test.o),
       cp = cp_ls, siap = out, pre_par = pre_par)
}
