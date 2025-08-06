siap.wrapper <- function(data, job, instance, id, ...) { 
  library(doParallel)
  library(Matrix)
  source("code/siap.R")
  x <- data
  
  d1 <- nrow(x)
  d2 <- ncol(x)
  
  S.valid.w <- which(instance$ind_valid.w == id) # may contain scattered missing entries in data
  # S.cal.w <- which(instance$ind_valid.w == id%%10 + 1)
  
  S.valid.o <- which(instance$ind_valid.o == id, arr.ind = T) # S.valid.o are all observed in data
  S.valid.o <- matrix(S.valid.o, ncol = 2)
  if (length(S.valid.w) > 0) {
    S.valid.o <- S.valid.o[-which(S.valid.o[,2] %in% S.valid.w), ] # 2: col
  }
  # S.valid.o2w <- S.valid.o[which(S.valid.o[,"col"] %in% S.cal.w), ]
  # S.valid.o <- S.valid.o[-which(S.valid.o[,"col"] %in% S.cal.w), ]
  
  x[S.valid.o] <- NA
  x[, S.valid.w] <- NA
  
  # S.cal.o <- which(instance$ind_valid.o == id%%10 + 1, arr.ind = T)
  # S.cal.o <- S.cal.o[-which(S.cal.o[,"col"] %in% S.cal.w), ]
  # S.cal.o <- S.cal.o[-which(S.cal.o[,"col"] %in% S.valid.w), ]
  
  M.tr <- !is.na(x)
  # M.tr[, S.cal.w] <- F
  M.tr[, S.valid.w] <- F
  # M.tr[S.cal.o] <- F
  M.tr[S.valid.o] <- F
  
  x.tr <- x
  x.tr[which(M.tr == F)] <- NA
  
  # -------------- fit model --------------
  
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
  
  # ---------------- cp ----------------
  # Ome.cal.o <- matrix(F, nrow(x), ncol(x))
  # Ome.cal.o[S.cal.o] <- T
  
  # cp_q.o <- function(alpha, i) {
  #   # S.cal.o_i <- S.cal.o[which(S.cal.o[,"row"]==i), , drop = F]
  #   # rsdl <- c(abs(x[S.cal.o_i] - x_imp[S.cal.o_i]))
  #   S.cal.o_i <- which(Ome.cal.o[i, ])
  #   rsdl <- c(abs(x[i, S.cal.o_i] - x_imp[i, S.cal.o_i]))
  #   return(quantile(rsdl, 1 - alpha, type = 1))
  # }
  
  # cp_q.w <- function(alpha, i) {
  #   rsdl <- c(abs(x[i, S.cal.w] - x_imp[i, S.cal.w]))
  #   return(quantile(rsdl, 1 - alpha, type = 1, na.rm = T))
  # }
  
   
  # setdiff(which(apply(x, 2, function(v) all(is.na(v)))), which(apply(data, 2, function(v) all(is.na(v))))) 
  
  # S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
  # S.test.o2w <- S.test.o[which(S.test.o[,"col"] %in% S.cal.w), ]
  # S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.cal.w), ]
  
  # Ome.valid <- is.na(x)&(!is.na(data))
  # Ome.valid.o <- matrix(F, nrow(x), ncol(x))
  # Ome.valid.o[S.valid.o] <- T
  # 
  # Ome.valid.w <- matrix(F, nrow(x), ncol(x)) # w and o2w pixels
  # Ome.valid.w[,S.valid.w] <- T
  # Ome.valid.w[S.valid.o2w] <- T
  
  # q.o <- sapply(1:nrow(x), function(i) cp_q.o(cp.alpha, i))
  # q.w <- sapply(1:nrow(x), function(i) cp_q.w(cp.alpha, i))
  
  # x_q <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  # for (i in 1:nrow(S.valid.o2w)) {
  #   x_q[S.valid.o2w[i,,drop=F]] <- q.w[S.valid.o2w[i,1]]
  # }
  # for (i in 1:nrow(S.valid.o)) {
  #   x_q[S.valid.o[i,,drop=F]] <- q.o[S.valid.o[i,1]]
  # }
  # x_q[,S.valid.w] <- q.w
  # 
  # ttl_coverage.w <- sum((abs(x_imp - data)<=x_q)*Ome.valid.w)
  # ttl_coverage.o <- sum((abs(x_imp - data)<=x_q)*Ome.valid.o)
  # ave_coverage_rate <- (ttl_coverage.w + ttl_coverage.o)/sum(Ome.valid)
  
  # cp_ls <- list(cp.alpha = cp.alpha, 
  #               S.cal.w = S.cal.w, S.cal.o = S.cal.o, S.valid.o2w = S.valid.o2w,
  #               cp_q.w=q.w, 
  #               cp_q.o=q.o)
  
  # ----------- rae-------------
  rr <- abs((x_imp - data)/data)
  
  #  calculated on test set but not S.test.o2w
  if (length(S.valid.w) == 0) {
    rae.w <- NA
  } else {
    rae.w <- rr[,S.valid.w,drop=F]
  }
  
  if (length(S.valid.o) == 0) {
    rae.o <- NA
  } else {
    mask.o <- matrix(NA, d1, d2)
    mask.o[S.valid.o] <- 1
    rae.o <- rr*mask.o
  }
  
  list(x = x, id = id, 
       mmrae.w = mean(rae.w, na.rm = T),
       mmrae.o = mean(rae.o, na.rm = T),
       # ave_coverage = ave_coverage_rate, 
       # ave_coverage.w = ttl_coverage.w/(sum(Ome.valid)-nrow(S.valid.o)), 
       # ave_coverage.o = ttl_coverage.o/nrow(S.valid.o),
       # cp = cp_ls, 
       siap = out)
}