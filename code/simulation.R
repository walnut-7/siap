library(batchtools)
source("code/read_data_functions.R")
source("code/preprocess_functions.R")
source("code/postprocess_functions.R")

source("code/simulation_siap_wrapper.R")
source("code/simulation_gp_wrapper.R")
source("code/simulation_marss_wrapper.R")


# ------- config -------
reg = try({
  loadRegistry(work.dir = getwd(), file.dir = "bt_simulation", writeable = T) # load registry
}, silent = T)
if (inherits(target, "try-error")) {
  reg = makeExperimentRegistry(work.dir = getwd(), file.dir = "bt_simulation", seed = 1) # create new registry
}

account = NA # change to slurm account name if you are running on a slurm cluster

# ------ data -------
dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt
d1 = nrow(xt)
d2 = ncol(xt)

# -------- helper functions ----------
{
  subsample = function(data, job, p.dt, p.sctt, cal = F, ...) {
    source("code/read_data_functions.R")
    d1 <- nrow(data)
    d2 <- ncol(data)
    mask <- chunk_mask(m = d1, n = d2, option = "uniform", ratio = p.dt, return = "mask")
    ind <- which(rbinom(d1*d2, 1, p.sctt) == 1)
    mask[ind] <- T
    return(which(mask == T)) # return missing entries
  }
  
  subsample_cal = function(data, job, p.dt, p.sctt, cal=F, p.cal, ...) {
    source("code/read_data_functions.R")
    d1 <- nrow(data)
    d2 <- ncol(data)
    mask <- chunk_mask(m = d1, n = d2, option = "uniform", ratio = p.dt, return = "mask")
    ind <- which(rbinom(d1*d2, 1, p.sctt) == 1)
    mask[ind] <- T
    
    # cp 
    if (cal == T) {
      S_full <- which(apply(mask, 2, function(v) all(v == T))) # fully missing
      # p.w <- length(S_full) / ncol(mask)
      # p.o <- (sum(mask == T) - length(S_full) * nrow(mask)) / ((ncol(mask) - length(S_full)) * nrow(mask))
      S.cal.w <- setdiff(1:ncol(mask), S_full)[as.logical(rbinom(ncol(mask) - length(S_full), 1, p.cal))]
      M.o <- (mask==F)
      rownames(M.o) <- 1:nrow(mask)
      M.o[, S.cal.w] <- F
      S.cal.o <- which(M.o == T, arr.ind = T)[as.logical(rbinom(sum(M.o), 1, p.cal)), ]
      return(list(miss = which(mask == T), cp = list(S.cal.o = S.cal.o, S.cal.w = S.cal.w))) 
    }
    return(which(mask == T)) # return missing entries
  }
  
  
  subsample_divide = function(data, job, p.dt, p.sctt, id, nv, ...) {
    source("code/read_data_functions.R")
    d1 <- nrow(data)
    d2 <- ncol(data)
    mask <- chunk_mask(m = d1, n = d2, option = "uniform", ratio = p.dt, return = "mask")
    ind <- which(rbinom(d1*d2, 1, p.sctt) == 1)
    mask[ind] <- T
    
    nb <- ceiling(nrow(data) / nv)
    gpidx <- rep(1:nb, each = nv)[1:nrow(data)]
    # data[gpidx == id, ]
    return(list(miss = which(mask == T), div = which(gpidx == id))) # return missing entries
  }
  
  get_jobdf <- function(ids, jobs_per_chunk = 10) {
    njobs <- nrow(ids)
    jobdf <- data.frame(job.id = ids$job.id, 
                        chunk=rep(1:ceiling(njobs / jobs_per_chunk), jobs_per_chunk)[1:njobs])
    return(jobdf)
  }
}

# ------- add problems and algorithms --------
addProblem(name = "ssi", data = xt, fun = subsample, reg = reg, seed = 1) 
addProblem(name = "ssi_cal", data = xt, fun = subsample_cal, reg = reg, seed = 1) 

addAlgorithm(name = "siap", fun = siap.wrapper, reg = reg)
addAlgorithm(name = "gp", fun = gp.wrapper, reg = reg)

addProblem(name = "ssi_batch", data = xt, fun = subsample_divide, reg = reg, seed = 1) 
addAlgorithm(name = "marss", fun = marss.wrapper, reg = reg)

# -------- submit jobs --------
{
  ## ----- add jobs: siap and gp -----
  siap_par <- data.frame(matrix(nrow = 7, ncol = 17))
  colnames(siap_par) <- c("cp.alpha", 
                          "step1", "step2", "step1.prd", "step2.ar", "step2.prd",
                          "rank.max", "step2.rank", "lambda", "lambda1", "lambda2", "alpha", "time.lag",
                          "use.ecm","cov.method", "thresh.sig", "rankSig")
  rownames(siap_par) <- c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "s1_si_trend")
  siap_par$p.cal = 0.1
  siap_par$cp.alpha = 0.05
  siap_par$step2.rank = 10
  siap_par$rank.max = 10
  siap_par$cov.method = "diag+lowrank"
  siap_par$prd.type = sapply(1:nrow(siap_par), function(...) list(c("l","y")))
  
  siap_par["siap",c("step1", "step2", "step1.prd", "step2.ar", "step2.prd",
                    "lambda", "lambda1", "lambda2", "alpha", "time.lag",
                    "use.ecm","cov.method", "thresh.sig", "rankSig")] <- 
    data.frame(step1 = T, step2 = T, step1.prd = T, step2.ar = T, step2.prd = T,
               lambda = 5, lambda1 = 3, lambda2 = 3, alpha = 3, time.lag = 2,
               use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100)
  
  siap_par["si",c("step1", "step2", "step1.prd", "lambda", "use.ecm")] <- 
    data.frame(step1 = T, step2 = F, step1.prd = F, lambda = 5, use.ecm = F)
  siap_par["si_trend",c("step1", "step2", "step1.prd", "lambda", "use.ecm")] <- 
    data.frame(step1 = T, step2 = F, step1.prd = T, lambda = 5, use.ecm = F)
  siap_par["si_trend_cov",c("step1", "step2", "step1.prd", "lambda", "use.ecm",
                            "cov.method", "thresh.sig", "rankSig")] <- 
    data.frame(step1 = T, step2 = F, step1.prd = T, lambda = 5, use.ecm = T,
               cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100)
  siap_par["sia",c("step1", "step2", "step2.prd", "step2.ar", "lambda1", "lambda2", "alpha", "time.lag")] <- 
    data.frame(step1 = F, step2 = T, step2.prd = F, step2.ar = T, lambda1 = 3, lambda2 = 3, alpha=3, time.lag=2)
  siap_par["sia_trend",c("step1", "step2", "step2.prd", "step2.ar", "lambda1", "lambda2", "alpha", "time.lag")] <- 
    data.frame(step1 = F, step2 = T, step2.prd = T, step2.ar = T, lambda1 = 3, lambda2 = 3, alpha=3, time.lag=2)
  
  siap_par["s1-si_trend", c("step1", "step2", "step1.prd", "step2.ar", "step2.prd",
                            "lambda", "lambda1", "lambda2", 
                            "use.ecm","cov.method", "thresh.sig", "rankSig")] <-
    data.frame(step1 = T, step2 = T, step1.prd = T, step2.ar = F, step2.prd = T,
               lambda = 5, lambda1 = 3, lambda2 = 3, 
               use.ecm = T, cov.method = "diag+lowrank", thresh.sig = 1e-2, rankSig = 100)
  
  ### ----- add gp jobs only -----
  ades = list(gp = data.frame(covfun = "exponential_isotropic"))
  pdes = list(ssi = data.frame(p.dt = c(0.1), p.sctt = c(0.1)))
  addExperiments(pdes, ades, repls = 100, reg = reg) # add jobs
  
  ### ----- add siap jobs only -----
  ades = list(siap = siap_par)
  pdes = list(ssi_cal = data.frame(p.dt = c(0.5), p.sctt = c(0.5), cal = T, p.cal = 0.1))
  addExperiments(pdes, ades, repls = 100, reg = reg) # add jobs
  
  ## ----- add jobs: marss -----
  nv <- 10 # number of variables per batch
  nb <- ceiling(nrow(xt) / nv)
  
  ### ----- add marss jobs only -----
  ades = list(marss = data.frame(model = "low-rank", p = 2, r = 1))
  pdes = list(ssi_batch = data.frame(p.dt = c(0.5), p.sctt = c(0.5), id = 1:nb, nv = nv))
  addExperiments(pdes, ades, repls = 100, reg = reg) # add jobs
  
  ## ----- add job tags -----
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==T & step2==T & step1.prd==T & step2.ar==T & step2.prd==T & use.ecm==T)), tags = rownames(siap_par)[1])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==T & step2==F & step1.prd==F & use.ecm==F)), tags = rownames(siap_par)[2])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==T & step2==F & step1.prd==T & use.ecm==F)), tags = rownames(siap_par)[3])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==T & step2==F & step1.prd==T & use.ecm==T)), tags = rownames(siap_par)[4])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==F & step2==T & step2.ar==T & step2.prd==F)), tags = rownames(siap_par)[5])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==F & step2==T & step2.ar==T & step2.prd==T)), tags = rownames(siap_par)[6])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap", algo.pars = (step1==T & step2==T & step2.ar==F & step2.prd==T)), tags = rownames(siap_par)[7])
  addJobTags(findExperiments(prob.name = "ssi_cal", algo.name = "siap"), tags = "current")
  
  addJobTags(findExperiments(prob.name = "ssi_batch", algo.name = "marss", algo.pars = (model == "low-rank")), tags = "low_rank")
  
  ## ----- submit jobs ------
  gp_ids <- getJobTable()$job.id[getJobTable()$problem == "ssi" & getJobTable()$algorithm %in% c("gp")]
  siap_ids <- ijoin(findExperiments(prob.name = "ssi_cal"), findTagged("current"))
  marss_ids <- findExperiments(prob.name = "ssi_batch", algo.pars = (model == "low-rank"))$job.id
  ids_ls <- list(gp = gp_ids, siap = siap_ids, marss = marss_ids)
  resources_ls <- list(gp = list(account = account, walltime = '8:00:00', memory='10000m', ncpus=8), # gp
                       siap = list(account = account, walltime = '8:00:00', memory='10000m', ncpus=10), # siap
                       marss = list(account = account, walltime = '12:00:00', memory='1000m', ncpus=4)) # MARSS
  for (i in 1:3) {
    submit_ids = ids_ls[[i]]
    # submit_ids = findErrors()
    resources = resources_ls[[i]] 
    resources$chunks.as.arrayjobs = TRUE # each unit job will have the resources specified by `resources`
    
    submitJobs(get_jobdf(submit_ids), resources = resources)
  }
}

# ------ wait for job completion -----
waitForJobs()
# getStatus()

# ------ concatenate marss blocks -------
target <- try({
  loadRegistry(work.dir = getwd(), file.dir = "bt_simu_marss", make.default = FALSE)
}, silent = T)
if (inherits(target, "try-error")) {
  target = makeRegistry(work.dir = getwd(), file.dir = "bt_simu_marss", seed = 1, make.default = FALSE)
  batchExport(export = list(xt = xt), reg = target)
  concat_marss <- function(repl, pdt) {
    library(batchtools)
    library(dplyr)
    reg = loadRegistry(work.dir = getwd(), file.dir = "bt_simulation", writeable = F)
    fn <- function(obj1, obj2) {
      x <- rbind(obj1$x, obj2$x)
      x_pred <- rbind(obj1$x_pred, obj2$x_pred)
      x_sd <- rbind(obj1$x_sd, obj2$x_sd)
      rae = list(rae.w=rbind(obj1$rae$rae.w, obj2$rae$rae.w), rae.o=rbind(obj1$rae$rae.o, obj2$rae$rae.o))
      return(list(x=x, x_pred=x_pred, x_sd=x_sd, rae=rae))
    }
    return(reduceResults(fn, findExperiments(prob.name = "ssi_batch", algo.name = "marss", prob.pars = (p.dt==pdt), repls = repl)$job.id))
  }
  
  batchMap(fun=concat_marss, args = data.table::CJ(repl = 1:100, pdt = c(0.1,0.3,0.5)), reg=target)
  resources = list(account = account, walltime = '1:00:00', memory='1000m', ncpus=4)
  resources$chunks.as.arrayjobs = TRUE 
  jobs_per_chunk = 20
  submitted_jobs <- findNotDone(reg = target)
  njobs <- nrow(submitted_jobs)
  jobdf <- data.frame(job.id = submitted_jobs, chunk=rep(1:ceiling(njobs / jobs_per_chunk), jobs_per_chunk)[1:njobs])
  submitJobs(reg = target, resources = resources, ids = jobdf)
  # getStatus(reg=target)
  # killJobs(reg=target)
}

# ------ wait for job completion -----
waitForJobs()
# getStatus()

# ------ generate intermediate results -------
{
  source("code/simulation_intermediate_result.R")
}

