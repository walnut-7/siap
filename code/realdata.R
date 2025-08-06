library(batchtools)
source("code/read_data_functions.R")
source("code/preprocess_functions.R")
source("code/postprocess_functions.R")
source("code/siap.R")

source("code/realdata_siap_wrapper.R")
source("code/realdata_gp_wrapper.R")
source("code/realdata_marss_wrapper.R")

# ------- config -------
reg = try({
  loadRegistry(work.dir = getwd(), file.dir = "bt_realdata", writeable = T) # load registry
}, silent = T)
if (inherits(target, "try-error")) {
  reg = makeExperimentRegistry(work.dir = getwd(), file.dir = "bt_realdata", seed = 1) # create new registry
}

account = NA # change to slurm account name if you are running on a slurm cluster

# ------ data -------
dl <- load_ssi(option = "true")
x <- dl$x
wvl <- as.numeric(rownames(x))
t <- ymd(colnames(x))
d1 <- nrow(x)
d2 <- ncol(x)

# -------- helper functions ----------
{
  # subsample_null = function(data, job, ...) {
  #   return(NULL) 
  # }
  
  subsample_cal <- function(data, job, p.valid, cal=F, p.cal, ...) {
    d1 <- nrow(data)
    d2 <- ncol(data)
    S_full <- which(apply(data, 2, function(x) all(is.na(x)))) # fully missing
    S.valid.w <- setdiff(1:d2, S_full)[as.logical(rbinom(d2 - length(S_full), 1, p.valid))]
    M.o <- !is.na(data) # is observed
    M.o[, S.valid.w] <- F
    S.valid.o <- which(M.o == T, arr.ind = T)[as.logical(rbinom(sum(M.o), 1, p.valid)), ]
    M.o[S.valid.o] <- F
    miss <- which(M.o==F)
    
    # ---------- cp -----------
    if (cal == T) {
      S_full <- which(apply(M.o, 2, function(v) all(v == F))) # fully missing
      S.cal.w <- setdiff(1:ncol(M.o), S_full)[as.logical(rbinom(ncol(M.o) - length(S_full), 1, p.cal))]
      M.o[, S.cal.w] <- F
      S.cal.o <- which(M.o == T, arr.ind = T)[as.logical(rbinom(sum(M.o), 1, p.cal)), ]
      return(list(miss = miss, cp = list(S.cal.o = S.cal.o, S.cal.w = S.cal.w))) 
    }
    
    return(list(miss = miss)) # return missing entries
  }
  
  subsample <- function(data, job, p.valid, ...) {
    d1 <- nrow(data)
    d2 <- ncol(data)
    
    S_full <- which(apply(data, 2, function(x) all(is.na(x)))) # fully missing
    
    S.valid.w <- setdiff(1:d2, S_full)[as.logical(rbinom(d2 - length(S_full), 1, p.valid))]
    M.o <- !is.na(data) # is observed
    M.o[, S.valid.w] <- F
    
    S.valid.o <- which(M.o == T, arr.ind = T)[as.logical(rbinom(sum(M.o), 1, p.valid)), ]
    
    M.o[S.valid.o] <- F
    
    return(which(M.o==F)) # return missing entries
  }
  
  subsample_divide = function(data, job, p.valid, id, nv, ...) {
    d1 <- nrow(data)
    d2 <- ncol(data)
    
    S_full <- which(apply(data, 2, function(x) all(is.na(x)))) # fully missing
    
    S.valid.w <- setdiff(1:d2, S_full)[as.logical(rbinom(d2 - length(S_full), 1, p.valid))]
    M.o <- !is.na(data)
    M.o[, S.valid.w] <- F
    
    S.valid.o <- which(M.o == T, arr.ind = T)[as.logical(rbinom(sum(M.o), 1, p.valid)), ]
    
    M.o[S.valid.o] <- F
    
    nb <- ceiling(nrow(data) / nv)
    gpidx <- rep(1:nb, each = nv)[1:nrow(data)]
    
    return(list(miss = which(M.o == F), div = which(gpidx == id))) # return missing entries
  }
  
  get_jobdf <- function(ids, jobs_per_chunk = 10) {
    njobs <- nrow(ids)
    jobdf <- data.frame(job.id = ids$job.id, 
                        chunk=rep(1:ceiling(njobs / jobs_per_chunk), jobs_per_chunk)[1:njobs])
    return(jobdf)
  }
}

# ------- add problems and algorithms --------
# addProblem(name = "ssi", data = x, fun = subsample_null, reg = reg, seed = 1) 
addProblem(name = "ssi_tr", data = x, fun = subsample, reg = reg, seed = 1) 
addProblem(name = "ssi_tr_batch", data = x, fun = subsample_divide, reg = reg, seed = 1) 
addProblem(name = "ssi_tr_cal", data = x, fun = subsample_cal, reg = reg, seed = 1) 

addAlgorithm(name = "gp", fun = gp.wrapper, reg = reg)
addAlgorithm(name = "siap", fun = siap.wrapper, reg = reg)
addAlgorithm(name = "marss", fun = marss.wrapper, reg = reg)

# -------- submit jobs --------
{
  ## --------- add jobs -----------
  ### --------- siap -----------
  cmd <- sprintf('Rscript hyper_param_tuning_siap.R %s', account) # hyper-parameter tuning
  system(cmd)
  
  lambda <- readRDS(file = "./output/realdata/lambda.rds")
  lambda1 <- readRDS(file = "./output/realdata/lambda1.rds")
  lambda2 <- readRDS(file = "./output/realdata/lambda1.rds")
  alpha <- readRDS(file = "./output/realdata/alpha.rds")
  p <- readRDS(file = "./output/realdata/p.rds")
  
  siap_par <- data.table::CJ(sorted = F, cp.alpha = 0.05,
                             step1 = T, step2 = T, step1.prd = T, step2.ar = T, step2.prd = T, 
                             prd.type = list(c("l","y"), c("l", "s", "y")),
                             rank.max = min(d1, d2), step2.rank = min(d1, d2), 
                             lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, time.lag = p,
                             use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100)
  ades = list(siap = siap_par)
  pdes = list(ssi_tr_cal = data.frame(p.valid = 0.1, cal = T, p.cal = 0.1))
  addExperiments(prob.designs = pdes, algo.designs = ades, repls = 1, reg = reg)
  
  
  ### -------------- gp -----------------
  ades = list(gp = data.frame(covfun = "exponential_isotropic"))
  pdes = list(ssi_tr = data.frame(p.valid = 0.1))
  addExperiments(prob.designs = pdes, algo.designs = ades, repls = 1, reg = reg)
  
  ### -------------- marss ---------------
  cmd <- sprintf('Rscript hyper_param_tuning_marss.R %s', account) # hyper-parameter tuning
  system(cmd)
  
  nv <- 10 # number of variables per batch
  nb <- ceiling(nrow(x) / nv)
  p <- readRDS(file = "./output/realdata/marss_p.rds")
  ades = list(marss = data.frame(model = "low-rank", p = p, r = 1))
  pdes = list(ssi_tr_batch = data.frame(p.valid = 0.1, id = 1:nb, nv = nv))
  
  ## ----- submit jobs ------
  resources = list(account = account, walltime = '4:00:00', memory='10000m', ncpus=10)
  submitJobs(findExperiments(prob.name = "ssi_tr_cal", 
                             algo.pars = (lambda == siap_par$lambda[2] & 
                                            lambda1 == siap_par$lambda1[2] & 
                                            lambda2 == siap_par$lambda2[2] & 
                                            alpha == siap_par$alpha[2] & 
                                            time.lag == siap_par$time.lag[2])) %>% ijoin(findNotDone()), resources = resources) # siap
  
  resources = list(account = account, walltime = '4:00:00', memory='10000m', ncpus=6)
  submitJobs(getJobTable()$job.id[getJobTable()$problem == "ssi_tr" & getJobTable()$algorithm %in% c("gp")], resources = resources) # gp
  
  resources = list(account = account, walltime = '4:00:00', memory='1000m', ncpus=10, chunks.as.arrayjobs = TRUE)
  submit_ids = findExperiments(prob.name = "ssi_tr_batch", algo.name = "marss", 
                               algo.pars = (p == ades$marss$p & r == ades$marss$r))
  submitJobs(get_jobdf(submit_ids, jobs_per_chunk = 4), resources = resources) # marss
}

# ------ wait for job completion -----
waitForJobs()
# getStatus()

# ------ generate intermediate results -------
{
  source("code/realdata_intermediate_result.R")
}

