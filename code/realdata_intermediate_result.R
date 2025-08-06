# ------------ prepare --------------
library(batchtools)
library(dplyr)
library(abind)
source("code/read_data_functions.R")
source("code/siap.R")

reg = try({
  loadRegistry(work.dir = getwd(), file.dir = "bt_realdata", writeable = T) # load registry
}, silent = T)
if (inherits(target, "try-error")) {
  reg = makeExperimentRegistry(work.dir = getwd(), file.dir = "bt_realdata", seed = 1) # create new registry
}

# ---------------------- siap and gp --------------------
c <- readRDS(file = "./output/realdata/lambda.rds")
a <- readRDS(file = "./output/realdata/lambda1.rds")
b <- readRDS(file = "./output/realdata/lambda1.rds")
d <- readRDS(file = "./output/realdata/alpha.rds")
p <- readRDS(file = "./output/realdata/p.rds")
ids <- findExperiments(prob.name = "ssi_tr_cal", 
                       algo.pars = (lambda == c & lambda1 == a & lambda2 == b & alpha == d & time.lag == p))
# siap0.res <- loadResult(
#   ids[which(sapply(unwrap(getJobPars(ids))$prd.type, function(v) identical(v, c("l","y")))), ]
# )
siap.res <- loadResult(
  ids[which(sapply(unwrap(getJobPars(ids))$prd.type, function(v) identical(v, c("l","s","y")))), ]
)

gp.res <- loadResult(
  findExperiments(prob.name = "ssi_tr", algo.name = "gp")
)

saveRDS(siap.res, file = "./output/realdata/siap_traintest.rds")
saveRDS(gp.res, file = "./output/realdata/gp_traintest.rds")

# -------------------- marss --------------------
custom_combine <- function(a, b) {
  x <- rbind(a$x, b$x)
  x_pred <- rbind(a$x_pred, b$x_pred)
  x_sd <- rbind(a$x_sd, b$x_sd)
  rae.w <- rbind(a$rae$rae.w , b$rae$rae.w)
  rae.o <- rbind(a$rae$rae.o , b$rae$rae.o)
  x_boots <- abind(a$x_boots, b$x_boots, along = 1)
  ttl_coverage.w <- sum(a$ttl_coverage.w, b$ttl_coverage.w)
  ttl_coverage.o <- sum(a$ttl_coverage.o, b$ttl_coverage.o)
  n.w <- sum(a$n.w, b$n.w)
  n.o <- sum(a$n.o, b$n.o)
  
  runtime.marss <- c(a$runtime$runtime.marss, b$runtime$runtime.marss)
  runtime.boot <- c(a$runtime$runtime.boot, b$runtime$runtime.boot)
  
  runtime = list(runtime.marss=runtime.marss, runtime.boot=runtime.boot)
  
  rae = list(rae.w=rae.w, rae.o=rae.o)
  
  list(x = x, x_pred = x_pred, x_sd = x_sd, rae = rae, x_boots = x_boots,
       ttl_coverage.w = ttl_coverage.w, n.w = n.w,
       ttl_coverage.o = ttl_coverage.o, n.o = n.o,
       runtime = runtime)
}

p <- readRDS(file = "./output/realdata/marss_p.rds")
ades = list(marss = data.frame(model = "low-rank", p = p, r = 1))
ids = findExperiments(prob.name = "ssi_tr_batch", algo.name = "marss", 
                      algo.pars = (p == ades$marss$p & r == ades$marss$r))
marss.res = reduceResults(fun=custom_combine, ids = ids, reg = reg)

saveRDS(marss.res, file = "./output/realdata/marss_lowrank_traintest.rds")
