compute_cp <- function(x, x_imp, xt, cp.alpha, S.cal.w, S.cal.o) {
  Ome.cal.o <- matrix(FALSE, nrow(x), ncol(x))
  Ome.cal.o[S.cal.o] <- TRUE

  cp_q.o <- function(alpha, i) {
    S.cal.o.i <- which(Ome.cal.o[i, ])
    rsdl <- c(abs(x[i, S.cal.o.i] - x_imp[i, S.cal.o.i]))
    quantile(rsdl, 1 - alpha, type = 1)
  }

  cp_q.w <- function(alpha, i) {
    rsdl <- c(abs(x[i, S.cal.w] - x_imp[i, S.cal.w]))
    quantile(rsdl, 1 - alpha, type = 1, na.rm = TRUE)
  }

  S.test.o <- which(is.na(x), arr.ind = TRUE)
  S.test.w <- which(apply(x, 2, function(v) all(is.na(v))))
  S.test.o <- S.test.o[-which(S.test.o[, "col"] %in% S.test.w), , drop = FALSE]
  S.test.o2w <- S.test.o[which(S.test.o[, "col"] %in% S.cal.w), , drop = FALSE]
  S.test.o <- S.test.o[-which(S.test.o[, "col"] %in% S.cal.w), , drop = FALSE]

  Ome.test <- is.na(x)
  Ome.test.o <- matrix(FALSE, nrow(x), ncol(x))
  Ome.test.o[S.test.o] <- TRUE

  Ome.test.w <- matrix(FALSE, nrow(x), ncol(x))
  Ome.test.w[, S.test.w] <- TRUE
  Ome.test.w[S.test.o2w] <- TRUE

  q.o <- sapply(1:nrow(x), function(i) cp_q.o(cp.alpha, i))
  q.w <- sapply(1:nrow(x), function(i) cp_q.w(cp.alpha, i))

  x_q <- matrix(0, nrow = nrow(xt), ncol = ncol(xt))
  for (i in 1:nrow(S.test.o2w)) {
    x_q[S.test.o2w[i, , drop = FALSE]] <- q.w[S.test.o2w[i, 1]]
  }
  for (i in 1:nrow(S.test.o)) {
    x_q[S.test.o[i, , drop = FALSE]] <- q.o[S.test.o[i, 1]]
  }
  x_q[, S.test.w] <- q.w

  covered <- abs(x_imp - xt) <= x_q
  ttl_coverage.w <- sum(covered * Ome.test.w)
  ttl_coverage.o <- sum(covered * Ome.test.o)
  ave_coverage_rate <- (ttl_coverage.w + ttl_coverage.o) / sum(is.na(x))

  wvl <- as.numeric(rownames(xt))
  coverage_wvl <- data.frame(
    wvl = wvl,
    overall = apply(covered * Ome.test, 1, function(v) sum(v)) / apply(Ome.test, 1, sum),
    o = apply(covered * Ome.test.o, 1, function(v) sum(v)) / apply(Ome.test.o, 1, sum),
    w = apply(covered * Ome.test.w, 1, function(v) sum(v)) / apply(Ome.test.w, 1, sum)
  )

  list(
    cp = list(
      cp.alpha = cp.alpha,
      S.cal.w = S.cal.w,
      S.cal.o = S.cal.o,
      S.test.o2w = S.test.o2w,
      cp_q.w = q.w,
      cp_q.o = q.o,
      coverage_wvl = coverage_wvl
    ),
    x_q = x_q,
    S.test.o = S.test.o,
    S.test.w = S.test.w,
    S.test.o2w = S.test.o2w,
    Ome.test = Ome.test,
    Ome.test.o = Ome.test.o,
    Ome.test.w = Ome.test.w,
    ttl_coverage.w = ttl_coverage.w,
    ttl_coverage.o = ttl_coverage.o,
    ave_coverage = ave_coverage_rate
  )
}
