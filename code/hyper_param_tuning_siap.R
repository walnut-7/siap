# Select hyper-parameters for siap.
args <- commandArgs(trailingOnly = TRUE)
account = args[1] # slurm account name

# ---------- prepare ---------
library(batchtools)
library(dplyr)
source("code/read_data_functions.R")
source("code/realdata_siap_wrapper_cv.R")

reg = try({
  loadRegistry(work.dir = getwd(), file.dir = "bt_realdata_cv", writeable = T) # load registry
}, silent = T)
if (inherits(target, "try-error")) {
  reg = makeExperimentRegistry(work.dir = getwd(), file.dir = "bt_realdata_cv", seed = 1) # create new registry
}

addAlgorithm(name = "siap", fun = siap.wrapper, reg = reg)

# ------------ data ----------
dl <- load_ssi(option = "true")
x <- dl$x
wvl <- as.numeric(rownames(x))
t <- ymd(colnames(x))
d1 <- nrow(x)
d2 <- ncol(x)

# ------------ helper functions ---------------
{
  subsample = function(data, job, n.fold, ...) {
    d1 <- nrow(data)
    d2 <- ncol(data)
    
    # cross validation sample splitting 
    S_full <- which(apply(data, 2, function(x) all(is.na(x)))) # fully missing
    ind_valid.w <- rep(0, d2) # 0: missing
    ind_valid.w[setdiff(1:d2, S_full)] <- sample(1:n.fold, size = d2 - length(S_full), replace = T)
    S_obs <- which(!is.na(data), arr.ind = T)
    ind_valid.o <- matrix(0, d1, d2) # 0: missing
    ind_valid.o[S_obs] <- sample(1:n.fold, size = nrow(S_obs), replace = T)
    
    return(list(ind_valid.o = ind_valid.o, ind_valid.w = ind_valid.w))
  }
  
  subsample_s1 = function(data, job, n.fold, ...) {
    d1 <- nrow(data)
    d2 <- ncol(data)
    
    # cross validation sample splitting 
    S_full <- which(apply(data, 2, function(x) all(is.na(x)))) # fully missing
    S_obs <- which(!is.na(data), arr.ind = T)
    ind_valid.o <- matrix(0, d1, d2) # 0: missing
    ind_valid.o[S_obs] <- sample(1:n.fold, size = nrow(S_obs), replace = T)
    
    return(list(ind_valid.o = ind_valid.o))
  }
  
  subsample_s2 = function(data, job, n.fold, ...) {
    d1 <- nrow(data)
    d2 <- ncol(data)
    
    # cross validation sample splitting 
    S_full <- which(apply(data, 2, function(x) all(is.na(x)))) # fully missing
    ind_valid.w <- rep(0, d2) # 0: missing
    ind_valid.w[setdiff(1:d2, S_full)] <- sample(1:n.fold, size = d2 - length(S_full), replace = T)
    
    return(list(ind_valid.w = ind_valid.w))
  }
  
  subsample_null <- function(data, job, ...) {
    return(list())
  }
  
  get_jobdf <- function(submit_ids, jobs_per_chunk = 10) {
    ids <- ajoin(submit_ids, findTagged("fixed_rank2"))
    njobs <- nrow(ids)
    jobdf <- data.frame(job.id = ids$job.id, 
                        chunk=rep(1:ceiling(njobs / jobs_per_chunk), jobs_per_chunk)[1:njobs])
    return(jobdf)
  }
}


# --------- 1: select S1 hyperparameter ------------
{
  ## --------- step 1 -------------
  addProblem(name = "s1_cv", data = x, fun = subsample_s1, reg = reg, seed = 1) 
  siap_par <- data.table::CJ(id = 1:10,
                             step1 = T, step2 = F, step1.prd = T, step2.ar = NA, step2.prd = NA, 
                             prd.type = list(c("l","y")),
                             rank.max = min(d1, d2), lambda = seq(1, 12.5, by = 0.1), lambda1 = NA, lambda2 = NA, alpha = NA, time.lag = NA,
                             use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100,
                             sorted = FALSE)
  ades = list(siap = siap_par)
  pdes = list(s1_cv = data.frame(n.fold = 10))
  addExperiments(pdes, ades, repls = 1, reg = reg)
  
  submit_ids = findExperiments(prob.name = "s1_cv") %>% ijoin(findNotDone())
  resources = list(account = account, walltime = '8:00:00', memory='10000m', ncpus=10) # siap

  resources$chunks.as.arrayjobs = TRUE 
  
  submitJobs(get_jobdf(submit_ids), resources = resources)
  
  waitForJobs()
  # getStatus()
  
  ## --------- s1 selection -------
  cat("Selecting S1 hyperparameter\n")
  cand <- (getJobPars(findExperiments(prob.name = "s1_cv")) %>% unwrap())[["lambda"]] %>% unique()
  err.valid <- rep(NA, length(cand))
  reduce_fn <- function(x, y) {
    mmrae.o <- x$mmrae.o + y$mmrae.o
    return(list(mmrae.o = mmrae.o))
  }
  for (i in 1:length(cand)) {
    a <- cand[i]
    ids <- findExperiments(prob.name = "s1_cv", algo.pars = (lambda == a))
    err.valid[i] <- reduceResults(reduce_fn, ids = ids)$mmrae.o / nrow(ids)
  }
  lambda <- cand[which.min(err.valid)]
  saveRDS(lambda, file = "./output/realdata/lambda.rds")
  saveRDS(data.frame(lambda = cand, err = err.valid), file = "./output/realdata/lambda_path.rds")
  
  ## --------- s1 retrain -------
  lambda <- readRDS(file = "./output/realdata/lambda.rds")
  addProblem(name = "s1_retrain", data = x, fun = subsample_null, reg = reg, seed = 1)
  siap_par <- data.table::CJ(id = 0,
                             step1 = T, step2 = F, step1.prd = T, step2.ar = NA, step2.prd = NA, 
                             prd.type = list(c("l","y")),
                             rank.max = min(d1, d2), lambda = lambda, lambda1 = NA, lambda2 = NA, alpha = NA, time.lag = NA,
                             use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100,
                             sorted = FALSE)
  ades = list(siap = siap_par)
  pdes = list(s1_retrain = data.frame())
  addExperiments(pdes, ades, repls = 1, reg = reg)
  resources = list(account = account, walltime = '1:00:00', memory='10000m', ncpus=10) 
  submitJobs(findExperiments(prob.name = "s1_retrain", algo.pars = (lambda == siap_par$lambda)), resources = resources)
  
  waitForJobs()
  # getStatus()
  
  res <- loadResult(findExperiments(prob.name = "s1_retrain", algo.pars = (lambda == siap_par$lambda)))
  x1 <- x
  S.test.o <- which(is.na(x), arr.ind = T)
  S.test.w <- which(apply(x, 2, function(v) all(is.na(v))))
  S.test.o <- S.test.o[-which(S.test.o[,2] %in% S.test.w), ]
  x1[S.test.o] <- res$siap$fit1$x_imp1[S.test.o]
  
}

# --------- 2: select S2.1 hyperparameter ------------
{
  ## --------- step 2.1 -------------
  addProblem(name = "s21_cv", data = x1, fun = subsample_s2, reg = reg, seed = 1) 
  lambda <- readRDS(file = "./output/realdata/lambda.rds")
  siap_par <- data.table::CJ(id = 1:10,
                             step1 = F, step2 = T, step1.prd = NA, step2.ar = F, step2.prd = T, 
                             prd.type = list(c("l","y")),
                             rank.max = min(d1, d2), step2.rank = min(d1, d2), lambda = lambda, lambda1 = seq(8,15, by = 0.1), lambda2 = NA, alpha = NA, time.lag = NA,
                             use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100,
                             sorted = FALSE)
  siap_par$lambda2 <- siap_par$lambda1
  
  ades = list(siap = siap_par)
  pdes = list(s21_cv = data.frame(n.fold = 10))
  addExperiments(pdes, ades, repls = 1, reg = reg)
  resources = list(account = account, walltime = '8:00:00', memory='10000m', ncpus=10, chunks.as.arrayjobs = T) # siap
  # submitJobs(get_jobdf(findExperiments(prob.name = "s21_cv")), resources = resources)
  submitJobs(get_jobdf(ijoin(findExperiments(prob.name = "s21_cv"), findNotStarted())), resources = resources)
  
  waitForJobs()
  # getStatus()
  
  ## --------- s2.1 selection -------
  cat("Selecting S2.1 hyperparameter\n")
  c <- readRDS(file = "./output/realdata/lambda.rds")
  cand <- (getJobPars(findExperiments(prob.name = "s21_cv", algo.pars = (id == 1)) %>% ajoin(findTagged("fixed_rank2"))) %>% unwrap() %>% filter(lambda == c))[, .(lambda1,lambda2)]
  err.valid <- rep(NA, nrow(cand))
  reduce_fn <- function(x, y) {
    mmrae.w <- x$mmrae.w + y$mmrae.w
    return(list(mmrae.w = mmrae.w))
  }
  for (i in 1:nrow(cand)) {
    a <- cand[i,1]
    b <- cand[i,2]
    ids <- findExperiments(prob.name = "s21_cv", algo.pars = (lambda1 == a & lambda2 == b & lambda == c))
    ids <- ajoin(ids, findTagged("fixed_rank2"))
    err.valid[i] <- reduceResults(reduce_fn, ids = ids)$mmrae.w / nrow(ids)
  }
  
  lambda <- cand[which.min(err.valid), ]
  saveRDS(lambda$lambda1, file = "./output/realdata/lambda1.rds")
  saveRDS(lambda$lambda2, file = "./output/realdata/lambda2.rds")
  saveRDS(data.frame(lambda1 = cand$lambda1, lambda2 = cand$lambda2, err = err.valid), file = "./output/realdata/lambda12_path.rds")
  
  ## --------- s21 retrain -------------
  lambda <- readRDS(file = "./output/realdata/lambda.rds")
  lambda1 <- readRDS(file = "./output/realdata/lambda1.rds")
  lambda2 <- readRDS(file = "./output/realdata/lambda2.rds")
  
  addProblem(name = "s21_retrain", data = x1, fun = subsample_null, reg = reg, seed = 1) 
  siap_par <- data.table::CJ(id = 0,
                             step1 = F, step2 = T, step1.prd = NA, step2.ar = F, step2.prd = T, 
                             prd.type = list(c("l","y")),
                             rank.max = min(d1, d2), step2.rank = min(d1, d2), lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, 
                             alpha = NA, time.lag = NA,
                             use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100,
                             sorted = FALSE)
  ades = list(siap = siap_par)
  pdes = list(s21_retrain = data.frame())
  addExperiments(pdes, ades, repls = 1, reg = reg)
  resources = list(account = account, walltime = '1:00:00', memory='10000m', ncpus=10) 
  submit_ids <- findExperiments(prob.name = "s21_retrain", algo.pars = (lambda1 == siap_par$lambda1 & lambda2 == siap_par$lambda2 & lambda == siap_par$lambda))
  if (nrow(submit_ids) == 0) {
    submit_ids <- findExperiments(prob.name = "s21_retrain", algo.pars = (lambda1 == siap_par$lambda1 & lambda2 == siap_par$lambda2)) # lambda = 10
  }
  submitJobs(submit_ids, resources = resources)
  
  waitForJobs()
  # getStatus()
  
}

# --------- 3: select p -------------
{
  a <- readRDS(file = "./output/realdata/lambda1.rds")
  b <- readRDS(file = "./output/realdata/lambda2.rds")
  c <- readRDS(file = "./output/realdata/lambda.rds")
  
  res <- loadResult(findExperiments(prob.name = "s21_retrain", algo.pars = (lambda1 == a & lambda2 == b & lambda == c)))
  tx = seq(1, d2)
  Phi = get_Phi(tx = tx,
                type = c("l","y"), spline.type = "b", 
                phase_l = lubridate::ymd(20191201) - lubridate::ymd(colnames(x)[1]),
                int.knots = list(l = c(778.50, 1605.25, 2445.50)),
                length.out = d2)
  b <- res$siap$fit2$b - Phi %*% res$siap$fit2$Theta 
  
  max_lag <- 12
  bic_vals <- numeric(max_lag)
  p_bic <- rep(NA, ncol(b))
  for (i in 1:ncol(b)) {
    y <- b[,i]
    for (p in 1:max_lag) {
      model <- arima(y, order = c(p, 0, 0), include.mean = TRUE, method = "ML")
      bic_vals[p] <- BIC(model)
    }
    p_bic[i] <- which.min(bic_vals)
  }
  # p <- floor(quantile(p_bic, 0.9))
  p <- max(p_bic)
  saveRDS(p, file = "./output/realdata/p.rds")
  saveRDS(p_bic, file = "./output/realdata/p_bic.rds")
}

# --------- 4: select S2.2 hyperparameter ------------
{
  ## --------- step 2.2 -------------
  lambda1 <- readRDS(file = "./output/realdata/lambda1.rds")
  lambda2 <- readRDS(file = "./output/realdata/lambda2.rds")
  lambda <- readRDS(file = "./output/realdata/lambda.rds")
  p <- readRDS(file = "./output/realdata/p.rds")
  addProblem(name = "s22_cv", data = x1, fun = subsample_s2, reg = reg, seed = 1) 

  siap_par <- data.table::CJ(id = 1:10,
                             step1 = F, step2 = T, step1.prd = NA, step2.ar = T, step2.prd = T, 
                             prd.type = list(c("l","y")),
                             rank.max = min(d1, d2), step2.rank = min(d1, d2), lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, 
                             alpha = seq(10, 3e6, by = 1e5), time.lag = p,
                             use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100,
                             sorted = FALSE)
  ades = list(siap = siap_par)
  pdes = list(s22_cv = data.frame(n.fold = 10))
  addExperiments(pdes, ades, repls = 1, reg = reg)
  resources = list(account = account, walltime = '8:00:00', memory='10000m', ncpus=10, chunks.as.arrayjobs = T) # siap
  submitJobs(get_jobdf(ijoin(findExperiments(prob.name = "s22_cv"), findNotDone())), resources = resources)
  
  waitForJobs()
  # getStatus()
  
  ## ----------- s2.2 selection -------------
  cat("Selecting S2.2 hyperparameter\n")
  a <- readRDS(file = "./output/realdata/lambda1.rds")
  b <- readRDS(file = "./output/realdata/lambda2.rds")
  c <- readRDS(file = "./output/realdata/lambda.rds")
  p <- readRDS(file = "./output/realdata/p.rds")
  cand <- (getJobPars(findExperiments(prob.name = "s22_cv") %>% ajoin(findTagged("fixed_rank2"))) %>% unwrap() %>% filter(lambda == c, lambda1 == a, lambda2 == b, time.lag == p))[["alpha"]] %>% unique()
  if (length(cand) == 0) {
    cand <- (getJobPars(findExperiments(prob.name = "s22_cv") %>% ajoin(findTagged("fixed_rank2"))) %>% unwrap() %>% filter(lambda1 == a, lambda2 == b, time.lag == p))[["alpha"]] %>% unique() # lambda = 10
  }
  err.valid <- rep(NA, length(cand))
  reduce_fn <- function(x, y) {
    mmrae.w <- x$mmrae.w + y$mmrae.w
    return(list(mmrae.w = mmrae.w))
  }
  for (i in 1:length(cand)) {
    d <- cand[i]
    ids <- findExperiments(prob.name = "s22_cv", algo.pars = (abs(alpha - d) < 1e-8 & lambda1 == a & lambda2 == b & lambda == c & time.lag == p))
    if (nrow(ids) == 0) {
      ids <- findExperiments(prob.name = "s22_cv", algo.pars = (abs(alpha - d) < 1e-8 & lambda1 == a & lambda2 == b & time.lag == p)) # lambda = 10
    }
    ids <- ids %>% ajoin(findTagged("fixed_rank2"))
    err.valid[i] <- reduceResults(reduce_fn, ids = ids)$mmrae.w / nrow(ids)
  }
  
  alpha <- cand[which.min(err.valid)]
  # plot(cand, err.valid)
  
  saveRDS(alpha, file = "./output/realdata/alpha.rds")
  saveRDS(data.frame(alpha = cand, err = err.valid), file = "./output/realdata/alpha_path.rds")
}

cat("SIAP hyper-parameter tuning completed.")