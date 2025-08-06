library(doParallel)

postprocess_func <- function(x_out, pre_par, x) {
  # Post-processing in a wrap. Consists of 3 steps:
  # 1. inverse normalization,
  # 2. inverse Box-Cox transformation, then minus b,
  # 3. cumulative sum.
  # Inputs:
  # - x_out: the data to be post-processed.
  # - pre_par: a list, the pre-processing parameters as output by preprocess_func().
  # - x: the original observed data matrix.
  
  # ---------- inverse normalization ----------
  x_out <- x_out * pre_par$nmlz_sd + pre_par$nmlz_mu
  # ---------- inverse Box-Cox transformation ------------
  if (pre_par$box.cox == F) {
    x_diff <- x_out
  } else {
    x_diff <- t(sapply(1:nrow(x_out), 
                       function(i) inv_bc(x_out[i,], pre_par$bc_lambda_ls[i], pre_par$bc_b_ls[i])))
  }
  # ---------- cumulative sum -------------
  x_imp <- t(sapply(1:nrow(x_diff), function(i) recover_from_diff(x_diff[i, ], x[i, ])))
  return(x_imp)
}

inv_bc <- function(y, lambda, b=0) {
  if (lambda == 0) {
    exp(y) - b
  } else if (lambda == 1) {
    y - b
  } else {
    (y * lambda + 1) ^ (1 / lambda) - b
  }
  
}

recover_from_diff <- function(v_diff, v, final.match = F) {
  # Recover the original data from its 1st order difference, 
  # given at least one 1 observation in the original vector.
  # - v_diff: the 1st order difference vector. Preferrably no NA's in v_diff.
  # - v: the original data vector, may contain NA's. length(v) == length(v_diff)+1
  anchor <- which(!is.na(v))[1]
  #anchor_and_after <- cumsum(c(v[anchor], v_diff[anchor:length(v_diff)]))
  if (anchor > 1) {
    anchor_and_before <- rev(cumsum(c(v[anchor], - rev(v_diff[1:(anchor - 1)]))))
  } else {
    anchor_and_before <- v[anchor]
  }
  #v_imp <- c(anchor_and_before, anchor_and_after[-1])
  
  # ------ alternative -------
  v_imp <- v
  if (anchor < length(v)) {
    for (i in (anchor+1):length(v)) {
      v_imp[i] <- v_imp[i-1] + v_diff[i-1]
    }
  }
  
  v_imp[1:anchor] <- anchor_and_before
  
  # final match with observed data
  if (final.match == T) {
    v_imp[!is.na(v)] = v[!is.na(v)]
  }
  
  return(v_imp)
}

