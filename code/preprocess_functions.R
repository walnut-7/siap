library(doParallel)

preprocess_func <- function(x, diff = T, box.cox = T, fudge = 0.2) {
  # Pre-processing in a wrap. Consists of 3 steps:
  # 1. (optional) take 1st order difference,
  # 2. (optional) plus b to make all values positive, then do Box-Cox transformation,
  # 3. normalize by rows.
  if (diff == T) {
    x_diff <- t(apply(x, 1, diff))
  } else {
    x_diff <- x
  }
  if (box.cox == T) {
    tt <- boxcox_transf(x_diff, fudge = 0.2)
    x_ready <- tt$x_diff_lam
    mu <- apply(x_ready, 1, mean, na.rm = T)
    sd <- apply(x_ready, 1, sd, na.rm = T)
    x_ready <- (x_ready - mu) / sd
    return(list(x_ready = x_ready, bc_lambda_ls = fudge_func(tt$lambda_ls, fudge = fudge), bc_b_ls = tt$b_ls,
                nmlz_mu = mu, nmlz_sd = sd, box.cox = box.cox))
  } else {
    x_ready <- x_diff
    mu <- apply(x_ready, 1, mean, na.rm = T)
    sd <- apply(x_ready, 1, sd, na.rm = T)
    x_ready <- (x_ready - mu) / sd
    return(list(x_ready = x_ready, nmlz_mu = mu, nmlz_sd = sd, diff = diff, box.cox = box.cox))
  }
  
  
}

fudge_func <- function(lambda, fudge = 0.2) {
  actual_lambda <- ifelse(abs(lambda) < fudge,
         0, 
         ifelse(abs(lambda-1) < fudge,
           1,
           lambda
           )
         )
  return(actual_lambda)
}

bc_func <- function(y, fudge = 0.2) {
  bc_trans <- caret::BoxCoxTrans(y, na.rm = T, fudge = 0.2)
  y_lam <- predict(bc_trans, y)
  lambda <- bc_trans$lambda
  #y_lam <- (y^lambda - 1) / lambda
  return(list(lambda = lambda, y_lam = y_lam))
}

boxcox_transf <- function(x, fudge = 0.2) {
  
  bc_func <- function(y, plotit = F) {
    bc_trans <- caret::BoxCoxTrans(y, na.rm = T, fudge = 0.2)
    y_lam <- predict(bc_trans, y)
    lambda <- bc_trans$lambda
    #y_lam <- (y^lambda - 1) / lambda
    return(list(lambda = lambda, y_lam = y_lam))
  }
  
  cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
  if(is.na(cores)) cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
  if(is.na(cores)) cores <- parallel::detectCores() - 1
  doParallel::registerDoParallel(cores)
  
  res <- foreach(j = 1:nrow(x), .combine="rbind", .multicombine=TRUE) %dopar% {
    y <- x[j,] 
    b <- - min(y, na.rm = T) + sd(y, na.rm = T)
    y <- y + b # ensure positive
    # y <- y - min(y, na.rm = T)
    # y <- y + sd(y, na.rm = T)
    bc <- bc_func(y = y)
    matrix(c(bc$lambda, b, bc$y_lam), nrow = 1)
  }
  
  lambda_ls <- res[,1]
  b_ls <- res[,2]
  x_diff_lam <- res[,-(1:2)]
  
  return(list(lambda_ls = lambda_ls, b_ls = b_ls, x_diff_lam = x_diff_lam))
}
