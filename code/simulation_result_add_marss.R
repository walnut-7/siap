# ---------- switch -----------
# 1: MRAE vs wvl and overall MRAE
# 2: interval length vs wvl
# 3: coverage vs wvl
# 4: overall coverage
# 5: model object
args <- commandArgs(trailingOnly = TRUE)
run_parts <- as.integer(strsplit(args[1], "\\s+")[[1]])

parse_bool_flag <- function(x, arg_name) {
  x_upper <- toupper(x)
  if (!x_upper %in% c("T", "F", "TRUE", "FALSE", "1", "0")) {
    stop(paste(arg_name, "must be one of T/F/TRUE/FALSE/1/0"))
  }
  x_upper %in% c("T", "TRUE", "1")
}

# beta.oneshot_flag <- TRUE
# if (length(args) >= 2) {
#   beta.oneshot_flag <- parse_bool_flag(args[2], "beta.oneshot")
# }

# siap_preprocess_mode <- "with"
# if (length(args) >= 3) {
#   siap_preprocess_mode <- tolower(args[3])
# }
# if (!siap_preprocess_mode %in% c("with", "without")) {
#   stop("siap preprocess mode must be one of: with, without")
# }

marss_diff_flag <- FALSE
if (length(args) >= 4) {
  marss_diff_flag <- parse_bool_flag(args[4], "marss diff")
}

gp_diff_flag <- FALSE
if (length(args) >= 5) {
  gp_diff_flag <- parse_bool_flag(args[5], "gp diff")
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

# ------------ prepare --------------
library(batchtools)
library(dplyr)
source("code/read_data_functions.R")
source("./code/registry_config.R")
source("./code/postprocess_functions.R")
dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt
rm(dl)
d1 = nrow(xt)
d2 = ncol(xt)
source_registry_location <- get_registry_location("bt_simulation")
marss_registry_location <- get_registry_location("bt_simu_marss")
reg = loadRegistry(
  work.dir = source_registry_location$work.dir,
  file.dir = source_registry_location$file.dir,
  writeable = F
)

# ----------- concatenate blocks -----------
target <- loadRegistry(
  work.dir = marss_registry_location$work.dir,
  file.dir = marss_registry_location$file.dir,
  make.default = FALSE
)
na_ids <- getJobPars(findErrors(reg=target),reg=target) %>% unwrap() # not done
done_ids <- getJobPars(findDone(reg=target),reg=target) %>% unwrap()

# ------------ helper functions ---------------
# marss_ids <- findExperiments(prob.name = "ssi_batch", algo.pars = (model == "low-rank"))$job.id
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

get_marss_ids <- function(p, diff_flag = marss_diff_flag) {
  ids <- findJobs(expr = (pdt==p), reg = target)
  filter_diff_tag(ids, reg = target, diff_flag = diff_flag)
}

# siap.ids <- ijoin(findExperiments(algo.name = "siap"), findTagged("current"))
# get_siap_ids <- function(algo, pdt, beta.oneshot = beta.oneshot_flag, preprocess_mode = siap_preprocess_mode) {
#   ids = findTagged(algo) %>% ijoin(siap.ids) %>% ijoin(findExperiments(prob.pars = (p.dt == pdt)))
#   if (beta.oneshot == F) {
#     ids <- ids %>% ijoin(findTagged("iterbeta"))
#   } else {
#     ids <- ids %>% anti_join(findTagged("iterbeta"), by = join_by(job.id))
#   }

#   if (algo %in% c("si", "sia")) {
#     ids <- ids %>% ijoin(findTagged("preprocess"))
#   } else if (preprocess_mode == "with") {
#     ids <- ids %>% ijoin(findTagged("preprocess"))
#   } else if (preprocess_mode == "without") {
#     ids <- ids %>% anti_join(findTagged("preprocess"), by = join_by(job.id))
#   }

#   return(ids)
# }

# -------------- 1: MRAE vs wvl and overall MRAE -----------
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

  for (pdt in c(0.1,0.3,0.5)) {
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
    
    mrae_stats <- readRDS(file = result_path("mrae_stats", pdt))
    mrae_stats <- rbind(mrae_stats, mrae_stats0)
    saveRDS(mrae_stats, file = result_path("mrae_stats", pdt))
    
    df1 <- res$mmrae
    df1$method = "marss"
    df <- readRDS(file = result_path("overall_rel_mrae_margin", pdt))
    df0 <- df %>% filter(method == 'siap') %>%
      mutate(mrae0 = mrae / (1 - rel_margin_mrae)) %>%
      select(type, mrae0)
    df1$rel_margin_mrae = (df0$mrae0 - df1$mrae) / df0$mrae0
    df <- rbind(df, df1)
    saveRDS(df, file = result_path("overall_rel_mrae_margin", pdt))
  }
  cat("Part 1 done\n")
}


# -------------- 2: interval length vs wvl -------------
if (2 %in% run_parts) {
  cat("Running Part 2: Interval length vs wvl\n")
  wvl <- as.numeric(rownames(xt))
  for (pdt in c(0.1,0.3,0.5)) {
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
    res <- reduceResultsList(fun=map_fn, ids = get_marss_ids(pdt), reg = target,
                             missing.val = data.frame(length.w = rep(NA, length(wvl)), length.o = rep(NA, length(wvl))))
    res <- do.call(rbind, res) %>% f1("marss")

    length_stats <- readRDS(file = result_path("length_stats", pdt))
    length_stats <- rbind(length_stats, res)
    saveRDS(length_stats, file = result_path("length_stats", pdt))
  }
  cat("Part 2 done\n")
}


# ------------- 3: coverage vs wvl ------------
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
    
    df <- readRDS(file = result_path("coverage_wvl_stats", pdt))
    df <- rbind(df, df1)
    saveRDS(df, file = result_path("coverage_wvl_stats", pdt))
    
  }
  cat("Part 3 done\n")
}

# ----------- 4: overall coverage --------
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
    ids <- get_marss_ids(pdt)
    res <- reduceResultsList(fun = map_fn, ids = ids, reg = target, 
                             missing.val = data.frame(overall = NA, o = NA, w = NA))
    res <-Reduce(reduce_fn, res)
    
    coverage_stats0 <- data.frame(mean_coverage = mean(res$overall, na.rm = T),
                                  sd_coverage  = sd(res$overall, na.rm = T),
                                  mean_coverage.w = mean(res$w, na.rm = T),
                                  sd_coverage.w = sd(res$w, na.rm = T),
                                  mean_coverage.o = mean(res$o, na.rm = T),
                                  sd_coverage.o = sd(res$o, na.rm = T),
                                  method = "marss")
    coverage_stats <- readRDS(file = result_path("coverage_stats", pdt))
    coverage_stats <- rbind(coverage_stats, coverage_stats0)
    saveRDS(coverage_stats, file = result_path("coverage_stats", pdt))
  }
  cat("Part 4 done\n")
}

# ------------ 5: model object -----------
if (5 %in% run_parts) {
  cat("Running Part 5: Model object\n")
  for (pdt in c(0.1,0.3,0.5)) {
    obj <- loadResult(id = get_marss_ids(pdt)[1], reg = target)
    saveRDS(obj, file = result_path("obj", pdt, "marss"))
  }
  cat("Part 5 done\n")
}
