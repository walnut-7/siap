siap.wrapper <- function(data, job, instance, cp.alpha, ...) { 
  library(doParallel)
  library(Matrix)
  source("code/siap.R")
  x <- data
  x[instance$miss] <- NA
  d1 <- nrow(x)
  d2 <- ncol(x)
  
  ## -------------- fit model --------------
  
  S.cal.w <- instance$cp$S.cal.w
  S.cal.o <- instance$cp$S.cal.o
  
  M.o <- !is.na(x)
  rownames(M.o) <- 1:nrow(x)
  M.o[, S.cal.w] <- F
  M.o[S.cal.o] <- F
  
  x.tr <- x
  x.tr[which(M.o == F)] <- NA
  
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
  
  # ----------- cp ------------  
  
  Ome.cal.o <- matrix(F, nrow(x), ncol(x))
  Ome.cal.o[S.cal.o] <- T
  
  cp_q.o <- function(alpha, i) {
    # S.cal.o_i <- S.cal.o[which(S.cal.o[,"row"]==i), , drop = F]
    # rsdl <- c(abs(x[S.cal.o_i] - x_imp[S.cal.o_i]))
    S.cal.o_i <- which(Ome.cal.o[i, ])
    rsdl <- c(abs(x[i, S.cal.o_i] - x_imp[i, S.cal.o_i]))
    return(quantile(rsdl, 1 - alpha, type = 1))
  }
  
  cp_q.w <- function(alpha, i) {
    rsdl <- c(abs(x[i, S.cal.w] - x_imp[i, S.cal.w]))
    return(quantile(rsdl, 1 - alpha, type = 1, na.rm = T))
  }
  
  S.test.o <- which(is.na(x), arr.ind = T)
  S.test.w <- which(apply(x, 2, function(x) all(is.na(x))))
  S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
  S.test.o2w <- S.test.o[which(S.test.o[,"col"] %in% S.cal.w), ]
  S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.cal.w), ]
  
  Ome.test <- is.na(x)
  Ome.test.o <- matrix(F, nrow(x), ncol(x))
  Ome.test.o[S.test.o] <- T
  
  Ome.test.w <- matrix(F, nrow(x), ncol(x)) # w and o2w pixels
  Ome.test.w[,S.test.w] <- T
  Ome.test.w[S.test.o2w] <- T
  
  q.o <- sapply(1:nrow(x), function(i) cp_q.o(cp.alpha, i))
  q.w <- sapply(1:nrow(x), function(i) cp_q.w(cp.alpha, i))
  
  x_q <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  for (i in 1:nrow(S.test.o2w)) {
    x_q[S.test.o2w[i,,drop=F]] <- q.w[S.test.o2w[i,1]]
  }
  for (i in 1:nrow(S.test.o)) {
    x_q[S.test.o[i,,drop=F]] <- q.o[S.test.o[i,1]]
  }
  x_q[,S.test.w] <- q.w
  
  ttl_coverage.w <- sum((abs(x_imp - data)<=x_q)*Ome.test.w) # w and o2w pixels
  ttl_coverage.o <- sum((abs(x_imp - data)<=x_q)*Ome.test.o)
  ave_coverage_rate <- (ttl_coverage.w + ttl_coverage.o)/sum(is.na(x))
  
  wvl <- as.numeric(rownames(data))
  df <- data.frame(wvl = wvl,
    overall = apply((abs(x_imp - data)<=x_q)*Ome.test, 1, function(v) sum(v)) / apply(Ome.test, 1, sum),
    o = apply((abs(x_imp - data)<=x_q)*Ome.test.o, 1, function(v) sum(v)) / apply(Ome.test.o, 1, sum),
    w = apply((abs(x_imp - data)<=x_q)*Ome.test.w, 1, function(v) sum(v)) / apply(Ome.test.w, 1, sum))
  
  cp_ls <- list(cp.alpha = cp.alpha, 
                S.cal.w = S.cal.w, S.cal.o = S.cal.o, S.test.o2w = S.test.o2w,
                cp_q.w=q.w, 
                cp_q.o=q.o,
                coverage_wvl = df)
  
  # ----------- spectral mrae-------------
  #  calculated on test set but not S.test.o2w
  rr <- (x_imp - data)/data
  mrae.w <- apply(rr[,S.test.w], 1, function(v) mean(abs(v)))
  
  mask.o <- matrix(NA, d1, d2)
  mask.o[S.test.o] <- 1
  mrae.o <- apply(rr*mask.o, 1, function(v) mean(abs(v), na.rm=T))
  
  list(x = x, mrae = list(mrae.w=mrae.w, mrae.o=mrae.o), 
       ave_coverage = ave_coverage_rate, 
       ave_coverage.w = ttl_coverage.w/(sum(is.na(x))-nrow(S.test.o)), 
       ave_coverage.o = ttl_coverage.o/nrow(S.test.o),
       cp = cp_ls, siap = out)
}