# Collect siap and gp results.
# Usage:
#   Rscript code/simulation_result.R [switches]


# ---------- switch -----------
# 1: MRAE vs wvl and overall MRAE
# 2: interval length vs wvl
# 3: coverage vs wvl
# 4: overall coverage
# 5: model object
args <- commandArgs(trailingOnly = TRUE)
run_parts <- as.integer(strsplit(args[3], "\\s+")[[1]])
add <- ifelse(args[1]==0, F, T)
model <- args[2]

parse_bool_flag <- function(x, arg_name) {
  x_upper <- toupper(x)
  if (!x_upper %in% c("T", "F", "TRUE", "FALSE", "1", "0")) {
    stop(paste(arg_name, "must be one of T/F/TRUE/FALSE/1/0"))
  }
  x_upper %in% c("T", "TRUE", "1")
}

beta.oneshot_flag <- T
if (length(args) >= 4) {
  beta.oneshot_flag <- parse_bool_flag(args[4], "beta.oneshot")
}

siap_preprocess_mode <- "with"
if (length(args) >= 5) {
  siap_preprocess_mode <- tolower(args[5])
}
if (!siap_preprocess_mode %in% c("with", "without")) {
  stop("siap preprocess mode must be one of: with, without")
}

gp_diff_flag <- FALSE
if (length(args) >= 6) {
  gp_diff_flag <- parse_bool_flag(args[6], "gp diff")
}

marss_diff_flag <- FALSE
if (length(args) >= 7) {
  marss_diff_flag <- parse_bool_flag(args[7], "marss diff")
}

diff_label <- function(flag) {
  if (flag) "diff" else "nodiff"
}

result_suffix <- paste0(
  "_gp_", diff_label(gp_diff_flag),
  "_marss_", diff_label(marss_diff_flag)
)

result_path <- function(stem, pdt, method = NULL) {
  file.path(
    "./output/simulation",
    paste0(
      if (is.null(method)) stem else paste0(method, "_", stem),
      result_suffix,
      "_", 
      pdt * 10,
      ".rds"
    )
  )
}
siap_methods <- c("siap", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "s1_si_trend")
all_methods <- c(siap_methods, "gp")
if (add==T && !model%in%all_methods) {
  stop(paste("Incorrect model specified:", model))
}

# ------------ prepare --------------
library(batchtools)
library(dplyr)
# library(doParallel)
# library(foreach)
source("./code/read_data_functions.R")
source("./code/registry_config.R")
source("./code/postprocess_functions.R")
invisible(gcinfo(verbose = FALSE))
dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt
rm(dl)
d1 = nrow(xt)
d2 = ncol(xt)
source_registry_location <- get_registry_location("bt_simulation")
gp_uq_registry_location <- get_registry_location("bt_simu_gp_uq")
reg = loadRegistry(
  work.dir = source_registry_location$work.dir,
  file.dir = source_registry_location$file.dir,
  writeable = F
)
# cl <- NULL
# n_cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
#   if(is.na(n_cores)) n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
#   doParallel::registerDoParallel(n_cores)
# cl <- parallel::makeCluster(n_cores)
# doParallel::registerDoParallel(cl)


# ------------ helper functions ---------------
retrieve_x <- function(obj, job) {
  x <- xt
  if (!is.atomic(job$instance)) { # siap
    x[job$instance$miss] <- NA
  } else { # gp
    x[job$instance] <- NA
  }
  return(x)
}

retrieve_imp <- function(obj, job, algo = c("siap", "gp")) {
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
        x = retrieve_x(obj, job)
        Ome.tr <- !is.na(x)
        Ome.tr[,obj$cp$S.cal.w] <- F
        Ome.tr[obj$cp$S.cal.o] <- F
        x_imp[Ome.tr] = x[Ome.tr]
      }
    }
  } else {
    x_imp = obj$x_pred
  }
  
  if (!is.null(obj$pre_par)) {
    x_imp <- postprocess_func(x_imp, obj$pre_par)
  }
  return(x_imp)
}

siap.ids <- ijoin(findExperiments(algo.name = "siap", reg = reg), findTagged("current", reg = reg))
get_siap_ids <- function(algo, pdt, beta.oneshot = beta.oneshot_flag, preprocess_mode = siap_preprocess_mode) {
  ids = findTagged(algo, reg = reg) %>% ijoin(siap.ids) %>% ijoin(findExperiments(prob.pars = (p.dt == pdt), reg = reg))
  if (beta.oneshot == F) {
    ids <- ids %>% ijoin(findTagged("iterbeta", reg = reg))
  } else {
    ids <- ids %>% anti_join(findTagged("iterbeta", reg = reg), by = join_by(job.id))
  }

  # "si" and "sia" are always evaluated with preprocess outputs.
  if (algo %in% c("si", "sia")) {
    ids <- ids %>% ijoin(findTagged("preprocess", reg = reg))
  } else if (preprocess_mode == "with") {
    ids <- ids %>% ijoin(findTagged("preprocess", reg = reg))
  } else if (preprocess_mode == "without") {
    ids <- ids %>% anti_join(findTagged("preprocess", reg = reg), by = join_by(job.id))
  }

  return(ids)
}

filter_diff_tag <- function(ids, reg, diff_flag) {
  ids_tbl <- if (is.data.frame(ids)) ids else data.frame(job.id = ids)
  diff_ids <- findTagged("diff", reg = reg)
  if (diff_flag) {
    ids_tbl <- ijoin(ids_tbl, diff_ids)
  } else {
    ids_tbl <- anti_join(ids_tbl, diff_ids, by = join_by(job.id))
  }
  ids_tbl$job.id
}

get_gp_ids <- function(pdt, reg, diff_flag = gp_diff_flag) {
  ids <- findExperiments(prob.name = "ssi", algo.name = "gp", prob.pars = (p.dt == pdt), reg = reg)
  filter_diff_tag(ids, reg = reg, diff_flag = diff_flag)
}

get_gp_uq_ids <- function(pdt, target, diff_flag = gp_diff_flag) {
  ids <- ijoin(findTagged("gp", reg = target), findTagged(as.character(pdt), reg = target))
  filter_diff_tag(ids, reg = target, diff_flag = diff_flag)
}

# -------------- 1: MRAE vs wvl and overall MRAE -----------

if (1 %in% run_parts) {
  cat("Running Part 1: MRAE vs wvl and overall MRAE\n")
  wvl <- as.numeric(rownames(xt))
  combine_func <- function(x, y) {
    mrae.w = rbind(x$mrae$mrae.w, y$mrae$mrae.w)
    mrae.o = rbind(x$mrae$mrae.o, y$mrae$mrae.o)    
    list(mrae = list(mrae.w = mrae.w, mrae.o = mrae.o))
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

  map_func1 <- function(res, job, algo = c("siap", "gp")) {
    algo = match.arg(algo)
    x <- retrieve_x(res, job)

    S.test.w <- which(apply(x, 2, function(v) all(is.na(v))))
    S.test.o <- which(is.na(x), arr.ind = T)
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    x_imp <- retrieve_imp(res, job, algo)
    rr = abs((xt - x_imp)/xt)
    mask.o <- matrix(NA, nrow(xt), ncol(xt))
    mask.o[S.test.o] <- 1
    
    mrae.w <- mean(rr[, S.test.w])
    mrae.o <- mean(rr*mask.o, na.rm = T)
    df <- data.frame(type = "w", mrae = mrae.w)
    df <- rbind(df, data.frame(type = "o", mrae = mrae.o))
    return(df)
  }
  map_func2 <- function(res, job, algo) {
    x <- retrieve_x(res, job)
    S.test.w <- which(apply(x, 2, function(v) all(is.na(v))))
    S.test.o <- which(is.na(x), arr.ind = T)
    S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
    
    x_imp <- t(res$siap$fit1$Theta) %*% t(res$siap$fit1$Phi)
    if (!is.null(res$pre_par)) {
      x_imp <- postprocess_func(x_imp, res$pre_par)
    }
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
    run_one_method <- function(method_name) {
      ids <- if (method_name == "gp") {
        get_gp_ids(pdt, reg = reg)
      } else {
        get_siap_ids(method_name, pdt)
      }
      reduceResults(fun = combine_func, ids = ids, reg = reg)
    }

    run_one_margin_method <- function(method_name) {
      ids <- if (method_name == "gp") {
        get_gp_ids(pdt, reg = reg)
      } else {
        get_siap_ids(method_name, pdt)
      }
      res <- if (method_name == "gp") {
        do.call("rbind", reduceResultsList(fun = map_func1, ids = ids, algo = "gp", reg = reg))
      } else {
        do.call("rbind", reduceResultsList(fun = map_func1, ids = ids, reg = reg))
      }
      res$method <- method_name
      res
    }

    if (add == F) {
      method <- all_methods
      # ------ wvl mrae -------
      # res_lists <- foreach::foreach(
      #   method_name = method,
      #   .packages = c("batchtools", "dplyr"),
      #   .combine = "c",
      #   .multicombine = TRUE
      # ) %do% {
      #   list(run_one_method(method_name) %>% f1(method = method_name))
      # }
      res_lists <- vector("list", length = length(method))
      for (i in seq_along(method)) {
        method_name <- method[i]
        res_lists[[i]] <- run_one_method(method_name) %>% f1(method = method_name)
      }
      
      gc()
      
      mrae_stats <- Reduce(rbind, res_lists)
      saveRDS(mrae_stats, file = result_path("mrae_stats", pdt))
      
      gc()
      # ------ overall mrae -------
      # margin_list <- foreach::foreach(
      #   method_name = method,
      #   .packages = c("batchtools", "dplyr"),
      #   .combine = "c",
      #   .multicombine = TRUE
      # ) %do% {
      #   list(run_one_margin_method(method_name))
      # }
      margin_list <- vector("list", length = length(method))
      for (i in seq_along(method)) {
        method_name <- method[i]
        margin_list[[i]] <- run_one_margin_method(method_name)
      }
      df <- Reduce(rbind, margin_list)
      
      df0 <- do.call("rbind", reduceResultsList(fun=map_func2, ids = get_siap_ids("si_trend",pdt), reg = reg))
      
      df$rel_margin_mrae = (df0$mrae - df$mrae) / df0$mrae
      saveRDS(df, file = result_path("overall_rel_mrae_margin", pdt))
    } else {
      # ------ wvl mrae -------
      mrae_stats <- readRDS(file = result_path("mrae_stats", pdt))
      if (model %in% mrae_stats$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- run_one_method(model)
      mrae_stats <- rbind(mrae_stats, f1(res, model))
      saveRDS(mrae_stats, file = result_path("mrae_stats", pdt))

      # ------ overall mrae -------
      df <- readRDS(file = result_path("overall_rel_mrae_margin", pdt))
      if (model %in% df$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- run_one_margin_method(model)
      res$method <- model
      df0 <- do.call("rbind", reduceResultsList(fun=map_func2, ids = get_siap_ids("si_trend",pdt), reg = reg))
      res$rel_margin_mrae = (df0$mrae - res$mrae) / df0$mrae
      df <- rbind(df, res)
      saveRDS(df, file = result_path("overall_rel_mrae_margin", pdt))
    }
  }
  cat("Part 1 done\n")
}

gc()

# -------------- 2: interval length vs wvl -------------
if (2 %in% run_parts) {
  cat("Running Part 2: Interval length vs wvl\n")
  wvl <- as.numeric(rownames(xt))
  combine_func <- function(res) {
    length.w = res$cp$cp_q.w
    length.o = res$cp$cp_q.o
    
    data.frame(length.w = length.w, length.o = length.o)
  }
  map_func2 <- function(res) {
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
    run_one_method <- function(method_name) {
      if (method_name == "gp") {
        res <- reduceResultsList(
          fun = map_func2,
          ids = get_gp_ids(pdt, reg = reg),
          reg = reg
        )
      } else {
        res <- reduceResultsList(fun = combine_func, ids = get_siap_ids(method_name, pdt), reg = reg)
      }
      do.call(rbind, res)
    }

    if (add == F) {
      method <- all_methods
      # res_lists <- foreach::foreach(
      #   method_name = method,
      #   .packages = c("batchtools", "dplyr"),
      #   .combine = "c",
      #   .multicombine = TRUE
      # ) %do% {
      #   list(run_one_method(method_name))
      # }
      res_lists <- vector("list", length = length(method))
      for (i in seq_along(method)) {
        method_name <- method[i]
        res_lists[[i]] <- run_one_method(method_name) %>% f1(method = method_name)
      }
      
      gc()
      
      length_stats <- Reduce(rbind, res_lists)
      saveRDS(length_stats, file = result_path("length_stats", pdt))
    } else {
      length_stats <- readRDS(file = result_path("length_stats", pdt))
      if (model %in% length_stats$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- run_one_method(model)
      
      length_stats <- rbind(length_stats, f1(res, model))
      saveRDS(length_stats, file = result_path("length_stats", pdt))
    }
  }
  cat("Part 2 done\n")
}

gc()

# ------------- 3: coverage vs wvl ------------
obj <- loadResult(get_siap_ids("si",0.1)[1], reg = reg)
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
  target = loadRegistry(
    work.dir = gp_uq_registry_location$work.dir,
    file.dir = gp_uq_registry_location$file.dir,
    writeable = F,
    make.default = FALSE
  )
  
  for (pdt in c(0.1,0.3,0.5)) {
    run_one_method <- function(method_name) {
      if (method_name == "gp") {
        res <- reduceResults(
          fun = f1,
          ids = get_gp_uq_ids(pdt, target = target),
          reg = target
        )$df
        temp_func(res, "gp")
      } else {
        res <- reduceResultsList(fun = get_coverage_vs_wvl, ids = get_siap_ids(method_name, pdt), reg = reg) %>% Reduce(f = reduce_fn)
        temp_func(res, method_name)
      }
    }

    if (add == F) {
      method <- all_methods
      # res_lists <- foreach::foreach(
      #   method_name = method,
      #   .packages = c("batchtools", "dplyr"),
      #   .combine = "c",
      #   .multicombine = TRUE
      # ) %do% {
      #   list(run_one_method(method_name))
      # }
      res_lists <- vector("list", length = length(method))
      for (i in seq_along(method)) {
        method_name <- method[i]
        res_lists[[i]] <- run_one_method(method_name)
      }
      df <- do.call(rbind, res_lists)
      saveRDS(df, file = result_path("coverage_wvl_stats", pdt))
    } else {
      df <- readRDS(file = result_path("coverage_wvl_stats", pdt))
      
      if (model %in% df$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- run_one_method(model)
      df <- rbind(df, res)
      saveRDS(df, file = result_path("coverage_wvl_stats", pdt))
    }
  }
  cat("Part 3 done\n")
}

gc()
# ---------- 4: overall coverage -----------
if (4 %in% run_parts) {
  cat("Running Part 4: Overall coverage\n")

  combine_func <- function(x, y) {    
    ave_coverage = c(x$ave_coverage, y$ave_coverage)
    ave_coverage.w = c(x$ave_coverage.w, y$ave_coverage.w)
    ave_coverage.o = c(x$ave_coverage.o, y$ave_coverage.o)
    
    list(ave_coverage = ave_coverage, ave_coverage.w = ave_coverage.w, ave_coverage.o = ave_coverage.o)
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
    run_one_method <- function(method_name) {
      ids <- if (method_name == "gp") {
        get_gp_ids(pdt, reg = reg)
      } else {
        get_siap_ids(method_name, pdt)
      }
      reduceResults(fun = combine_func, ids = ids, reg = reg)
    }

    if (add == F) {
      method <- all_methods
      run_one_method <- function(method_name) {
        ids <- if (method_name == "gp") {
          get_gp_ids(pdt, reg = reg)
        } else {
          get_siap_ids(method_name, pdt)
        }
        reduceResults(fun = combine_func, ids = ids, reg = reg)
      }
      # res_lists <- foreach::foreach(
      #   method_name = method,
      #   .packages = c("batchtools", "dplyr"),
      #   .combine = "c",
      #   .multicombine = TRUE
      # ) %do% {
      #   list(run_one_method(method_name))
      # }
      res_lists <- vector("list", length = length(method))
      for (i in seq_along(method)) {
        method_name <- method[i]
        res_lists[[i]] <- run_one_method(method_name)
      }
      
      gc()
            
      coverage_stats <- do.call(rbind, lapply(seq_along(method), function(i) f2(res = res_lists[[i]], method = method[i])))
      saveRDS(coverage_stats, file = result_path("coverage_stats", pdt))
    } else {
      coverage_stats <- readRDS(file = result_path("coverage_stats", pdt))
      if (model %in% coverage_stats$method) {
        stop(paste("Model", model, "result already exist."))
      }
      res <- run_one_method(model)
      coverage_stats <- rbind(coverage_stats, f2(res, model))
      saveRDS(coverage_stats, file = result_path("coverage_stats", pdt))
    }
  }
  cat("Part 4 done\n")
}

gc()

# ------------ 5: model object -----------

if (5 %in% run_parts) {
  cat("Running Part 5: Model object\n")
  for (pdt in c(0.1,0.3,0.5)) {
    if (add == F) {
      method <- all_methods
      
      siap.obj <- loadResult(id = get_siap_ids("siap",pdt)[1], reg = reg)
      saveRDS(siap.obj, file = result_path("obj", pdt, method[1]))
      
      si.obj <- loadResult(id = get_siap_ids("si", pdt)[1], reg = reg)
      saveRDS(si.obj, file = result_path("obj", pdt, method[2]))
      
      si_trend.obj <- loadResult(id = get_siap_ids("si_trend",pdt)[1], reg = reg)
      saveRDS(si_trend.obj, file = result_path("obj", pdt, method[3]))
      
      si_trend_cov.obj <- loadResult(id = get_siap_ids("si_trend_cov", pdt)[1], reg = reg)
      saveRDS(si_trend_cov.obj, file = result_path("obj", pdt, method[4]))
      
      sia.obj <- loadResult(id = get_siap_ids("sia",pdt)[1], reg = reg)
      saveRDS(sia.obj, file = result_path("obj", pdt, method[5]))
      
      sia_trend.obj <- loadResult(id = get_siap_ids("sia_trend",pdt)[1], reg = reg)
      saveRDS(sia_trend.obj, file = result_path("obj", pdt, method[6]))
      
      s1_si_trend.obj <- loadResult(id = get_siap_ids("s1_si_trend",pdt)[1], reg = reg)
      saveRDS(s1_si_trend.obj, file = result_path("obj", pdt, method[7]))
      
      gp.obj <- loadResult(id = get_gp_ids(pdt, reg = reg)[1], reg = reg)
      saveRDS(gp.obj, file = result_path("obj", pdt, method[8]))
    } else {
      if (model == "gp") {
        obj <- loadResult(id = get_gp_ids(pdt, reg = reg)[1], reg = reg)
      } else {
        obj <- loadResult(id = get_siap_ids(model,pdt)[1], reg = reg)
      }
      saveRDS(obj, file = result_path("obj", pdt, model))
    }
  }
  cat("Part 5 done\n")
}

# if (!is.null(cl)) {
#   parallel::stopCluster(cl)
# }
