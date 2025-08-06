# --------- siap and gp results -----------

## ---------- switch -----------
# 1: MRAE vs wvl and overall coverage
# 2: interval length vs wvl
# 3: coverage vs wvl
# 4: overall mrae
# 5: model object

# args <- commandArgs(trailingOnly = TRUE)
# run_parts <- as.integer(strsplit(args[3], "\\s+")[[1]])

run_parts <- c(1, 2, 3, 4, 5)

add <- F # T: add one model result to the current result
model <- NULL # name of the model being added
if (add==T && !model%in%c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "s1_si_trend")) {
  stop(paste("Incorrect model specified:", model))
} else if (add==T) {
  cat(paste("Adding", model, "model Part", paste0(run_parts, collapse = ", "),"results."))
} else {
  cat(paste("Generating", "Part", paste0(run_parts, collapse = ", "),"results for siap and gp models."))
}

## ------------ prepare --------------
library(batchtools)
library(dplyr)
source("code/read_data_functions.R")

dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt
d1 = nrow(xt)
d2 = ncol(xt)
reg = loadRegistry(work.dir = getwd(), file.dir = "bt_simulation", writeable = F)

## ------------ helper functions ---------------
retrieve_imp <- function(obj, algo = c("siap", "gp")) {
  algo = match.arg(algo)
  if (algo == "siap") {
    out = obj$siap
    x_imp = out$fit2$x_imp
    if (is.null(x_imp)) {
      x_imp = out$fit1$x_imp1
      if (is.null(x_imp)) {
        if (out$fit1$prd == T) {
          x_imp = out$fit1$a %*% t(out$fit1$b)+ t(out$fit1$Theta) %*% t(out$fit1$Phi)
        } else {
          x_imp = out$fit1$a %*% t(out$fit1$b)
        }
        Ome.tr <- !is.na(obj$x)
        Ome.tr[,obj$cp$S.cal.w] <- F
        Ome.tr[obj$cp$S.cal.o] <- F
        x_imp[Ome.tr] = obj$x[Ome.tr]
      }
    }
  } else {
    x_imp = obj$x_pred
  }
  return(x_imp)
}

siap.ids <- ijoin(findExperiments(algo.name = "siap"), findTagged("current"))
get_siap_ids <- function(algo, pdt) {
  ids = findTagged(algo) %>% ijoin(siap.ids) %>% ijoin(findExperiments(prob.pars = (p.dt == pdt)))
  return(ids)
}

## -------------- 1: MRAE vs wvl and overall coverage -----------

if (1 %in% run_parts) {
  cat("Running Part 1: MRAE vs wvl and overall coverage\n")
  wvl <- as.numeric(rownames(xt))
  combine_func <- function(x, y) {
    mrae.w = rbind(x$mrae$mrae.w, y$mrae$mrae.w)
    mrae.o = rbind(x$mrae$mrae.o, y$mrae$mrae.o)
    
    ave_coverage = c(x$ave_coverage, y$ave_coverage)
    ave_coverage.w = c(x$ave_coverage.w, y$ave_coverage.w)
    ave_coverage.o = c(x$ave_coverage.o, y$ave_coverage.o)
    
    list(mrae = list(mrae.w = mrae.w, mrae.o = mrae.o), 
         ave_coverage = ave_coverage, ave_coverage.w = ave_coverage.w, ave_coverage.o = ave_coverage.o)
  }
  f1 <- function(res, method) {
    mrae_w <- res$mrae$mrae.w
    mrae_o <- res$mrae$mrae.o
    
    mean_mrae_w <- apply(mrae_w, 2, mean)
    sd_mrae_w <- apply(mrae_w, 2, sd)
    
    mean_mrae_o <- apply(mrae_o, 2, mean)
    sd_mrae_o <- apply(mrae_o, 2, sd)
    
    data.frame(
      wvl = wvl,
      mean_mrae.w = mean_mrae_w,
      sd_mrae.w = sd_mrae_w,
      mean_mrae.o = mean_mrae_o,
      sd_mrae.o = sd_mrae_o,
      method = method
    )
  }
  f2 <- function(res, method) {
    data.frame(
      mean_coverage = mean(res$ave_coverage),
      sd_coverage = sd(res$ave_coverage),
      mean_coverage.w = mean(res$ave_coverage.w),
      sd_coverage.w = sd(res$ave_coverage.w),
      mean_coverage.o = mean(res$ave_coverage.o),
      sd_coverage.o = sd(res$ave_coverage.o),
      method = method
    )
  }
  
  for (pdt in c(0.1,0.3,0.5)) {
    if (add == F) {
      siap.res <- reduceResults(fun=combine_func, ids = get_siap_ids("siap", pdt))
      si.res <- reduceResults(fun=combine_func, ids = get_siap_ids("si", pdt))
      si_trend.res <- reduceResults(fun=combine_func, ids = get_siap_ids("si_trend", pdt))
      si_trend_cov.res <- reduceResults(fun=combine_func, ids = get_siap_ids("si_trend_cov",pdt))
      sia.res <- reduceResults(fun=combine_func, ids = get_siap_ids("sia",pdt))
      sia_trend.res <- reduceResults(fun=combine_func, ids = get_siap_ids("sia_trend",pdt))
      gp.res <- reduceResults(fun=combine_func, ids = findExperiments(prob.name = "ssi", algo.name = "gp", prob.pars = (p.dt == pdt)))
      
      gc()
      
      method <- c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "gp")
      res_lists <- list(siap.res, si.res, si_trend.res, si_trend_cov.res, sia.res, sia_trend.res, gp.res)
      
      mrae_stats <- do.call(rbind, lapply(1:7, function(i) f1(res = res_lists[[i]], method = method[i])))
      saveRDS(mrae_stats, file = paste0("./output/simulation/mrae_stats_",pdt*10,".rds"))
      
      coverage_stats <- do.call(rbind, lapply(1:7, function(i) f2(res = res_lists[[i]], method = method[i])))
      saveRDS(coverage_stats, file = paste0("./output/simulation/coverage_stats_",pdt*10,".rds"))
    } else {
      mrae_stats <- readRDS(file = paste0("./output/simulation/mrae_stats_",pdt*10,".rds"))
      if (model %in% mrae_stats$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- reduceResults(fun=combine_func, ids = get_siap_ids(model, pdt))
      mrae_stats <- rbind(mrae_stats, f1(res, model))
      saveRDS(mrae_stats, file = paste0("./output/simulation/mrae_stats_",pdt*10,".rds"))
      
      coverage_stats <- readRDS(file = paste0("./output/simulation/coverage_stats_",pdt*10,".rds"))
      if (model %in% coverage_stats$method) {
        stop(paste("Model", model, "result already exist."))
      }
      coverage_stats <- rbind(coverage_stats, f2(res, model))
      saveRDS(coverage_stats, file = paste0("./output/simulation/coverage_stats_",pdt*10,".rds"))
    }
  }
  cat("Part 1 done\n")
}


## -------------- 2: interval length vs wvl -------------
if (2 %in% run_parts) {
  cat("Running Part 2: Interval length vs wvl\n")
  wvl <- as.numeric(rownames(xt))
  combine_func <- function(res) {
    length.w = res$cp$cp_q.w
    length.o = res$cp$cp_q.o
    
    data.frame(length.w = length.w, length.o = length.o)
  }
  combine_func2 <- function(res) {
    x_sd = res$x_sd
    x_sd[!is.na(res$x)] = NA
    S.test.o <- which(is.na(res$x), arr.ind = T)
    S.test.w <- which(apply(res$x, 2, function(v) all(is.na(v))))
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    mean_length.w = apply(qnorm(0.975)*x_sd[, S.test.w, drop=F], 1, mean)
    x_sd[, S.test.w] = NA
    mean_length.o = apply(x_sd, 1, mean, na.rm = T)
    
    data.frame(length.w = mean_length.w, length.o = mean_length.o)
  }
  f1 <- function(res, method) {
    data.frame(
      wvl = wvl,
      mean_length.w = apply(matrix(res$length.w, nrow=length(wvl)), 1, mean),
      mean_length.o = apply(matrix(res$length.o, nrow=length(wvl)), 1, mean),
      method = method
    )
  }
  for (pdt in c(0.1,0.3,0.5)) {
    if (add == F) {
      siap.res <- reduceResultsList(fun=combine_func, ids = get_siap_ids("siap",pdt))
      siap.res <- do.call(rbind, siap.res)
      si.res <- reduceResultsList(fun=combine_func, ids = get_siap_ids("si",pdt))
      si.res <- do.call(rbind, si.res)
      si_trend.res <- reduceResultsList(fun=combine_func, ids = get_siap_ids("si_trend", pdt))
      si_trend.res <- do.call(rbind, si_trend.res)
      si_trend_cov.res <- reduceResultsList(fun=combine_func, ids = get_siap_ids("si_trend_cov", pdt))
      si_trend_cov.res <- do.call(rbind, si_trend_cov.res)
      sia.res <- reduceResultsList(fun=combine_func, ids = get_siap_ids("sia", pdt))
      sia.res <- do.call(rbind, sia.res)
      sia_trend.res <- reduceResultsList(fun=combine_func, ids = get_siap_ids("sia_trend", pdt))
      sia_trend.res <- do.call(rbind, sia_trend.res)
      
      gp.res <- reduceResultsList(fun=combine_func2, ids = findExperiments(prob.name = "ssi", algo.name = "gp", prob.pars = (p.dt == pdt)))
      gp.res <- do.call(rbind, gp.res)
      
      gc()
      
      method <- c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "gp")
      res_lists <- list(siap.res, si.res, si_trend.res, si_trend_cov.res, sia.res, sia_trend.res, gp.res)
      
      length_stats <- do.call(rbind, lapply(1:7, function(i) f1(res = res_lists[[i]], method = method[i])))
      saveRDS(length_stats, file = paste0("./output/simulation/length_stats_",pdt*10,".rds"))
    } else {
      length_stats <- readRDS(file = paste0("./output/simulation/length_stats_",pdt*10,".rds"))
      if (model %in% length_stats$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- reduceResultsList(fun=combine_func, ids = get_siap_ids(model, pdt))
      res <- do.call(rbind, res)
      
      length_stats <- rbind(length_stats, f1(res, model))
      saveRDS(length_stats, file = paste0("./output/simulation/length_stats_",pdt*10,".rds"))
    }
  }
  cat("Part 2 done\n")
}

## ------------- 3: coverage vs wvl ------------
obj <- loadResult(get_siap_ids("si",0.1)[1])
if (is.null(obj$cp$coverage_wvl)) {
  version = "old"
} else {
  version = "new"
}

if (3 %in% run_parts) {
  wvl <- as.numeric(rownames(xt))
  if (version == "old") {
    cat("Running Part 3 (old version): coverage vs wvl\n")
    get_coverage_vs_wvl <- function(obj){
      x_imp <- retrieve_imp(obj, algo = "siap")
      S.test.o <- which(is.na(obj$x), arr.ind = T)
      S.test.w <- which(apply(obj$x, 2, function(v) all(is.na(v))))
      S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
      S.test.o2w <- obj$cp$S.test.o2w # S.test.o[which(S.test.o[,"col"] %in% obj$cp$S.cal.w), ]
      S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% obj$cp$S.cal.w), ]
      
      q.o <- obj$cp$cp_q.o
      q.w <- obj$cp$cp_q.w
      
      x_q <- matrix(0, nrow = nrow(obj$x), ncol = ncol(obj$x))
      for (i in 1:nrow(S.test.o2w)) {
        x_q[S.test.o2w[i,,drop=F]] <- q.w[S.test.o2w[i,1]]
      }
      for (i in 1:nrow(S.test.o)) {
        x_q[S.test.o[i,,drop=F]] <- q.o[S.test.o[i,1]]
      }
      x_q[,S.test.w] <- q.w
      
      Ome.test <- is.na(obj$x)
      Ome.test.o <- matrix(F, nrow(obj$x), ncol(obj$x))
      Ome.test.o[S.test.o] <- T
      
      Ome.test.w <- matrix(F, nrow(obj$x), ncol(obj$x)) # w and o2w pixels
      Ome.test.w[,S.test.w] <- T
      Ome.test.w[S.test.o2w] <- T
      
      wvl <- as.numeric(rownames(xt))
      df <- data.frame(wvl = wvl,
                       overall = apply((abs(x_imp - xt)<=x_q)*Ome.test, 1, function(v) sum(v)) / apply(Ome.test, 1, sum),
                       o = apply((abs(x_imp - xt)<=x_q)*Ome.test.o, 1, function(v) sum(v)) / apply(Ome.test.o, 1, sum),
                       w = apply((abs(x_imp - xt)<=x_q)*Ome.test.w, 1, function(v) sum(v)) / apply(Ome.test.w, 1, sum))
      
      return(df)
    }
  } else {
    cat("Running Part 3 (new version): coverage vs wvl\n")
    get_coverage_vs_wvl <- function(res) {
      return(res$cp$coverage_wvl)
    }
  }
  
  reduce_fn <- function(x, y) {
    overall <- rbind(x$overall, y$overall)
    o <- rbind(x$o, y$o)
    w <- rbind(x$w, y$w)
    list(overall = overall, o = o, w = w)
  }
  
  temp_func <- function(res, method) {
    data.frame(wvl = wvl,
               overall = apply(res$overall, 2, mean, na.rm = T),
               sd_overall = apply(res$overall, 2, sd, na.rm = T),
               o = apply(res$o, 2, mean, na.rm = T),
               sd_o = apply(res$o, 2, sd, na.rm = T),
               w = apply(res$w, 2, mean, na.rm = T),
               sd_w = apply(res$w, 2, sd, na.rm = T),
               method = method)
  }
  
  f1 <- function(x, y) {
    overall <- rbind(x$df$overall, y$df$overall)
    o <- rbind(x$df$o, y$df$o)
    w <- rbind(x$df$w, y$df$w)
    df = list(overall = overall, o = o, w = w)
    return(list(df = df))
  }
  target = loadRegistry(work.dir = "/home/yuxuank/solar", file.dir = "/home/yuxuank/bt_simu_gp_uq", writeable = F, make.default = FALSE)
  
  for (pdt in c(0.1,0.3,0.5)) {
    if (add == F) {
      siap.res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids("siap",pdt)) %>% Reduce(f=reduce_fn) %>% temp_func("siap")
      si.res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids("si",pdt)) %>% Reduce(f=reduce_fn) %>% temp_func("si")
      si_trend.res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids("si_trend", pdt)) %>% Reduce(f=reduce_fn) %>% temp_func("si_trend")
      si_trend_cov.res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids("si_trend_cov", pdt)) %>% Reduce(f=reduce_fn) %>% temp_func("si_trend_cov")
      sia.res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids("sia", pdt)) %>% Reduce(f=reduce_fn) %>% temp_func("sia")
      sia_trend.res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids("sia_trend", pdt)) %>% Reduce(f=reduce_fn) %>% temp_func("sia_trend")
      gp.res <- reduceResults(fun = f1, ids = ijoin(findTagged("gp", reg = target), findTagged(as.character(pdt), reg = target)), reg = target)$df %>% temp_func("gp")
      gc()
      df <- rbind(siap.res, si.res, si_trend.res, si_trend_cov.res, sia.res, sia_trend.res, gp.res)
      saveRDS(df, file = paste0("./output/simulation/coverage_wvl_stats_",pdt*10,".rds"))
    } else {
      df <- readRDS(file = paste0("./output/simulation/coverage_wvl_stats_",pdt*10,".rds"))
      
      if (model %in% df$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- reduceResultsList(fun=get_coverage_vs_wvl, ids = get_siap_ids(model,pdt)) %>% Reduce(f=reduce_fn) %>% temp_func(model)
      df <- rbind(df, res)
      saveRDS(df, file = paste0("./output/simulation/coverage_wvl_stats_",pdt*10,".rds"))
    }
  }
  cat("Part 3 done\n")
}


## ---------- 4: overall mrae -----------
if (4 %in% run_parts) {
  cat("Running Part 4: Overall MRAE\n")
  combine_func <- function(res, algo = c("siap", "gp")) {
    algo = match.arg(algo)
    S.test.w <- which(apply(res$x, 2, function(v) all(is.na(v))))
    S.test.o <- which(is.na(res$x), arr.ind = T)
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    x_imp <- retrieve_imp(res, algo)
    rr = abs((xt - x_imp)/xt)
    mask.o <- matrix(NA, nrow(xt), ncol(xt))
    mask.o[S.test.o] <- 1
    
    mrae.w <- mean(rr[, S.test.w])
    mrae.o <- mean(rr*mask.o, na.rm = T)
    df <- data.frame(type = "w", mrae = mrae.w)
    df <- rbind(df, data.frame(type = "o", mrae = mrae.o))
    
    return(df)
  }
  combine_func2 <- function(res, algo) {
    S.test.w <- which(apply(res$x, 2, function(v) all(is.na(v))))
    S.test.o <- which(is.na(res$x), arr.ind = T)
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    x_imp <- t(res$siap$fit1$Theta) %*% t(res$siap$fit1$Phi)
    rr = abs((xt - x_imp)/xt)
    mask.o <- matrix(NA, nrow(xt), ncol(xt))
    mask.o[S.test.o] <- 1
    
    mrae.w <- mean(rr[, S.test.w])
    mrae.o <- mean(rr*mask.o, na.rm = T)
    df <- data.frame(type = "w", mrae = mrae.w)
    df <- rbind(df, data.frame(type = "o", mrae = mrae.o))
    
    return(df)
  }
  for (pdt in c(0.1,0.3,0.5)) {
    if (add == F) {
      method <- c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "gp")
      
      siap.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids("siap",pdt)))
      siap.res$method <- method[1]
      si.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids("si", pdt)))
      si.res$method <- method[2]
      si_trend.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids("si_trend",pdt)))
      si_trend.res$method <- method[3]
      si_trend_cov.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids("si_trend_cov", pdt)))
      si_trend_cov.res$method <- method[4]
      sia.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids("sia", pdt)))
      sia.res$method <- method[5]
      sia_trend.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids("sia_trend", pdt)))
      sia_trend.res$method <- method[6]
      gp.res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = findExperiments(prob.name = "ssi", algo.name = "gp", prob.pars = (p.dt == pdt)), algo = "gp"))
      gp.res$method <- method[7]
      df <- rbind(siap.res, si.res, si_trend.res, si_trend_cov.res, sia.res, sia_trend.res, gp.res)
      
      df0 <- do.call("rbind", reduceResultsList(fun=combine_func2, ids = get_siap_ids("si_trend",pdt)))
      
      df$rel_margin_mrae = (df0$mrae - df$mrae) / df0$mrae
      saveRDS(df, file = paste0("./output/simulation/overall_rel_mrae_margin_", 10*pdt, ".rds"))
    } else {
      df <- readRDS(file = paste0("./output/simulation/overall_rel_mrae_margin_", 10*pdt, ".rds"))
      
      if (model %in% df$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- do.call("rbind", reduceResultsList(fun=combine_func, ids = get_siap_ids(model,pdt)))
      res$method <- model
      
      df0 <- do.call("rbind", reduceResultsList(fun=combine_func2, ids = get_siap_ids("si_trend",pdt)))
      res$rel_margin_mrae = (df0$mrae - res$mrae) / df0$mrae
      
      df <- rbind(df, res)
      saveRDS(df, file = paste0("./output/simulation/overall_rel_mrae_margin_", 10*pdt, ".rds"))
    }
  }
  cat("Part 4 done\n")
}

## ------------ 5: model object -----------

if (5 %in% run_parts) {
  cat("Running Part 5: Model object\n")
  for (pdt in c(0.1,0.3,0.5)) {
    if (add == F) {
      method <- c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "gp")
      
      siap.obj <- loadResult(id = get_siap_ids("siap",pdt)[1])
      saveRDS(siap.obj, file = paste0("./output/simulation/",method[1],"_obj_",pdt*10,".rds"))
      
      si.obj <- loadResult(id = get_siap_ids("si", pdt)[1])
      saveRDS(si.obj, file = paste0("./output/simulation/",method[2],"_obj_",pdt*10,".rds"))
      
      si_trend.obj <- loadResult(id = get_siap_ids("si_trend",pdt)[1])
      saveRDS(si_trend.obj, file = paste0("./output/simulation/",method[3],"_obj_",pdt*10,".rds"))
      
      si_trend_cov.obj <- loadResult(id = get_siap_ids("si_trend_cov", pdt)[1])
      saveRDS(si_trend_cov.obj, file = paste0("./output/simulation/",method[4],"_obj_",pdt*10,".rds"))
      
      sia.obj <- loadResult(id = get_siap_ids("sia",pdt)[1])
      saveRDS(sia.obj, file = paste0("./output/simulation/",method[5],"_obj_",pdt*10,".rds"))
      
      sia_trend.obj <- loadResult(id = get_siap_ids("sia_trend",pdt)[1])
      saveRDS(sia_trend.obj, file = paste0("./output/simulation/",method[6],"_obj_",pdt*10,".rds"))
      
      gp.obj <- loadResult(id = findExperiments(prob.name = "ssi", algo.name = "gp", prob.pars = (p.dt == pdt))[1])
      saveRDS(gp.obj, file = paste0("./output/simulation/",method[7],"_obj_",pdt*10,".rds"))
    } else {
      obj <- loadResult(id = get_siap_ids(model,pdt)[1])
      saveRDS(obj, file = paste0("./output/simulation/", model ,"_obj_",pdt*10,".rds"))
    }
  }
  cat("Part 5 done\n")
}

# --------------- add marss result --------------

## --------------- switch ----------------
# 1: MRAE vs wvl and overall MRAE
# 2: interval length vs wvl
# 3: coverage vs wvl
# 4: overall coverage
# 5: model object
run_parts <- c(1, 2, 3, 4, 5)

## ------------- prepare -------------
target <- loadRegistry(work.dir = "/home/yuxuank/solar", file.dir = "/home/yuxuank/bt_simu_marss", make.default = FALSE)
na_ids <- getJobPars(findErrors(reg=target),reg=target) %>% unwrap() # not done
done_ids <- getJobPars(findDone(reg=target),reg=target) %>% unwrap()
cat(paste("Adding MARSS model Part", paste0(run_parts, collapse = ", "), "results."))

## ------------ helper functions ---------------
# marss_ids <- findExperiments(prob.name = "ssi_batch", algo.pars = (model == "low-rank"))$job.id
get_marss_ids <- function(p) {
  ids <- findJobs(expr = (pdt==p), reg = target)
  return(ids$job.id)
}
siap.ids <- ijoin(findExperiments(algo.name = "siap"), findTagged("current"))
get_siap_ids <- function(algo, pdt) {
  ids = findTagged(algo) %>% ijoin(siap.ids) %>% ijoin(findExperiments(prob.pars = (p.dt == pdt)))
  return(ids)
}

## -------------- 1: MRAE vs wvl and overall MRAE -----------
if (1 %in% run_parts) {
  cat("Running Part 1: MRAE vs wvl and overall MRAE\n")
  wvl <- as.numeric(rownames(xt))
  map_fn <- function(obj) {
    mrae.w = apply(obj$rae$rae.w, 1, mean)
    mrae.o = apply(obj$rae$rae.o, 1, mean, na.rm = T)
    
    mmrae.w = mean(obj$rae$rae.w)
    mmrae.o = mean(obj$rae$rae.o, na.rm = T)
    mmrae <- data.frame(type = "w", mrae = mmrae.w)
    mmrae <- rbind(mmrae, data.frame(type = "o", mrae = mmrae.o))
    
    
    list(mrae.w = mrae.w, mrae.o = mrae.o, 
         mmrae = mmrae)
  }
  
  reduce_fn <- function(x, y) {
    mrae.w <- rbind(x$mrae.w, y$mrae.w)
    mrae.o <- rbind(x$mrae.o, y$mrae.o)
    
    mmrae <- rbind(x$mmrae, y$mmrae)
    
    list(mrae.w = mrae.w, mrae.o = mrae.o, 
         mmrae = mmrae)
  }
  
  combine_func2 <- function(res, algo) {
    S.test.w <- which(apply(res$x, 2, function(v) all(is.na(v))))
    S.test.o <- which(is.na(res$x), arr.ind = T)
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    x_imp <- t(res$siap$fit1$Theta) %*% t(res$siap$fit1$Phi)
    rr = abs((xt - x_imp)/xt)
    mask.o <- matrix(NA, nrow(xt), ncol(xt))
    mask.o[S.test.o] <- 1
    
    mrae.w <- mean(rr[, S.test.w])
    mrae.o <- mean(rr*mask.o, na.rm = T)
    df <- data.frame(type = "w", mrae = mrae.w)
    df <- rbind(df, data.frame(type = "o", mrae = mrae.o))
    
    return(df)
  }
  
  
  for (pdt in c(0.1,0.3,0.5)) {
    mrae_stats <- readRDS(file = paste0("./output/simulation/mrae_stats_",pdt*10,".rds"))
    if ("marss" %in% mrae_stats$method) {
      stop(paste("Model", model, "result already exist."))
    }
    ids <- get_marss_ids(pdt)
    res <- reduceResultsList(fun = map_fn, ids = ids, reg = target, 
                             missing.val = list(mrae.w = NA, mrae.o = NA, mmrae = data.frame(type = c("w", "o"), mrae = NA)))
    res <-Reduce(reduce_fn, res)
    mrae_stats0 <- data.frame(wvl = wvl,
                              mean_mrae.w = apply(res$mrae.w, 2, mean, na.rm = T),
                              sd_mrae.w = apply(res$mrae.w, 2, sd, na.rm = T),
                              mean_mrae.o = apply(res$mrae.o, 2, mean, na.rm = T),
                              sd_mrae.o = apply(res$mrae.o, 2, sd, na.rm = T),
                              method = "marss")
    mrae_stats <- rbind(mrae_stats, mrae_stats0)
    saveRDS(mrae_stats, file = paste0("./output/simulation/mrae_stats_",pdt*10,".rds"))
    
    df <- readRDS(file = paste0("./output/simulation/overall_rel_mrae_margin_",pdt*10,".rds"))
    if ("marss" %in% df$method) {
      stop(paste("Model", model, "result already exist."))
    }
    df1 <- res$mmrae
    df1$method = "marss"
    df0 <- do.call("rbind", reduceResultsList(fun=combine_func2, ids = get_siap_ids("si_trend",pdt)))
    df1$rel_margin_mrae = (df0$mrae - df1$mrae) / df0$mrae
    df <- rbind(df, df1)
    saveRDS(df, file = paste0("./output/simulation/overall_rel_mrae_margin_",pdt*10,".rds"))
  }
  cat("Part 1 done\n")
}


## -------------- 2: interval length vs wvl -------------
if (2 %in% run_parts) {
  cat("Running Part 2: Interval length vs wvl\n")
  wvl <- as.numeric(rownames(xt))
  map_fn <- function(obj) {
    x_sd = obj$x_sd
    x_sd[!is.na(obj$x)] = NA
    S.test.o <- which(is.na(obj$x), arr.ind = T)
    S.test.w <- which(apply(obj$x, 2, function(v) all(is.na(v))))
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    mean_length.w = apply(qnorm(0.975)*x_sd[, S.test.w, drop=F], 1, mean)
    x_sd[, S.test.w] = NA
    mean_length.o = apply(x_sd, 1, mean, na.rm = T)
    
    data.frame(length.w = mean_length.w, length.o = mean_length.o)
  }
  f1 <- function(res, method) {
    data.frame(
      wvl = wvl,
      mean_length.w = apply(matrix(res$length.w, nrow=length(wvl)), 1, mean, na.rm = T),
      mean_length.o = apply(matrix(res$length.o, nrow=length(wvl)), 1, mean, na.rm = T),
      method = method
    )
  }
  for (pdt in c(0.1,0.3,0.5)) {
    length_stats <- readRDS(file = paste0("./output/simulation/length_stats_",pdt*10,".rds"))
    if ("marss" %in% length_stats$method) {
      stop(paste("Model", model, "result already exist."))
    }
    res <- reduceResultsList(fun=map_fn, ids = get_marss_ids(pdt), reg = target,
                             missing.val = data.frame(length.w = rep(NA, length(wvl)), length.o = rep(NA, length(wvl))))
    res <- do.call(rbind, res) %>% f1("marss")
    length_stats <- rbind(length_stats, res)
    saveRDS(length_stats, file = paste0("./output/simulation/length_stats_",pdt*10,".rds"))
  }
  cat("Part 2 done\n")
}


## ------------- 3: coverage vs wvl ------------
if (3 %in% run_parts) {
  cat("Running Part 3: Coverage vs wvl\n")
  wvl <- as.numeric(rownames(xt))
  map_fn <- function(obj){
    x_imp <- obj$x_pred
    S.test.o <- which(is.na(obj$x), arr.ind = T)
    S.test.w <- which(apply(obj$x, 2, function(v) all(is.na(v))))
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    x_q <- qnorm(0.975)*obj$x_sd
    
    Ome.test <- is.na(obj$x)
    Ome.test.o <- matrix(F, nrow(obj$x), ncol(obj$x))
    Ome.test.o[S.test.o] <- T
    
    Ome.test.w <- matrix(F, nrow(obj$x), ncol(obj$x)) 
    Ome.test.w[,S.test.w] <- T
    
    df <- data.frame(#wvl = wvl,
      overall = apply((abs(x_imp - xt)<=x_q)*Ome.test, 1, function(v) sum(v)) / apply(Ome.test, 1, sum),
      o = apply((abs(x_imp - xt)<=x_q)*Ome.test.o, 1, function(v) sum(v)) / apply(Ome.test.o, 1, sum),
      w = apply((abs(x_imp - xt)<=x_q)*Ome.test.w, 1, function(v) sum(v)) / apply(Ome.test.w, 1, sum))
    
    return(df)
  }
  
  reduce_fn <- function(x, y) {
    overall <- rbind(x$overall, y$overall)
    o <- rbind(x$o, y$o)
    w <- rbind(x$w, y$w)
    
    list(overall = overall, o = o, w = w)
  }
  
  
  for (pdt in c(0.1,0.3,0.5)) {
    df <- readRDS(file = paste0("./output/simulation/coverage_wvl_stats_",pdt*10,".rds"))
    if ("marss" %in% df$method) {
      stop(paste("Model", model, "result already exist."))
    }
    ids <- get_marss_ids(pdt)
    res <- reduceResultsList(fun = map_fn, ids = ids, reg = target, 
                             missing.val = data.frame(overall = NA, o = NA, w = NA))
    res <-Reduce(reduce_fn, res)
    
    df1 <- data.frame(wvl = wvl,
                      overall = apply(res$overall, 2, mean, na.rm = T),
                      sd_overall = apply(res$overall, 2, sd, na.rm = T),
                      o = apply(res$o, 2, mean, na.rm = T),
                      sd_o = apply(res$o, 2, sd, na.rm = T),
                      w = apply(res$w, 2, mean, na.rm = T),
                      sd_w = apply(res$w, 2, sd, na.rm = T),
                      method = "marss")
    df <- rbind(df, df1)
    saveRDS(df, file = paste0("./output/simulation/coverage_wvl_stats_",pdt*10,".rds"))
    
  }
  cat("Part 3 done\n")
}

## ----------- 4: overall coverage --------
if (4 %in% run_parts) {
  cat("Running Part 4: Overall coverage\n")
  wvl <- as.numeric(rownames(xt))
  map_fn <- function(obj){
    x_imp <- obj$x_pred
    S.test.o <- which(is.na(obj$x), arr.ind = T)
    S.test.w <- which(apply(obj$x, 2, function(v) all(is.na(v))))
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    x_q <- qnorm(0.975)*obj$x_sd
    
    Ome.test <- is.na(obj$x)
    Ome.test.o <- matrix(F, nrow(obj$x), ncol(obj$x))
    Ome.test.o[S.test.o] <- T
    
    Ome.test.w <- matrix(F, nrow(obj$x), ncol(obj$x)) 
    Ome.test.w[,S.test.w] <- T
    
    df <- data.frame(overall = sum((abs(x_imp - xt)<=x_q)*Ome.test) / sum(Ome.test),
                     o = sum((abs(x_imp - xt)<=x_q)*Ome.test.o) / sum(Ome.test.o),
                     w = sum((abs(x_imp - xt)<=x_q)*Ome.test.w) / sum(Ome.test.w))
    
    return(df)
  }
  
  reduce_fn <- function(x, y) {
    overall <- c(x$overall, y$overall)
    o <- c(x$o, y$o)
    w <- c(x$w, y$w)
    list(overall = overall, o = o, w = w)
  }
  
  for (pdt in c(0.1,0.3,0.5)){
    coverage_stats <- readRDS(file = paste0("./output/simulation/coverage_stats_",pdt*10,".rds"))
    if ("marss" %in% coverage_stats$method) {
      stop(paste("Model", model, "result already exist."))
    }
    ids <- get_marss_ids(pdt)
    res <- reduceResultsList(fun = map_fn, ids = ids, reg = target, 
                             missing.val = data.frame(overall = NA, o = NA, w = NA))
    res <-Reduce(reduce_fn, res)
    
    coverage_stats0 <- data.frame(mean_coverage = mean(res$overall, na.rm = T),
                                  sd_coverage  = sd(res$overall, na.rm = T),
                                  mean_coverage.w = mean(res$w, na.rm = T),
                                  sd_coverage.w = sd(res$w, na.rm = T),
                                  mean_coverage.o = mean(res$w, na.rm = T),
                                  sd_coverage.o = sd(res$w, na.rm = T),
                                  method = "marss")
    coverage_stats <- rbind(coverage_stats, coverage_stats0)
    saveRDS(coverage_stats, file = paste0("./output/simulation/coverage_stats_",pdt*10,".rds"))
  }
  cat("Part 4 done\n")
}

## ------------ 5: model object -----------
if (5 %in% run_parts) {
  cat("Running Part 5: Model object\n")
  for (pdt in c(0.1,0.3,0.5)) {
    obj <- loadResult(id = get_marss_ids(pdt)[1], reg = target)
    saveRDS(obj, file = paste0("./output/simulation/marss_obj_",pdt*10,".rds"))
  }
  cat("Part 5 done\n")
}


