# Usage:
#   Rscript code/simulation_postexp_gp_uq.R [target_file_dir]
#   target_file_dir <- "/path/to/bt_simu_gp_uq"; source("code/simulation_postexp_gp_uq.R")

library(dplyr)
library(batchtools)

source("./code/read_data_functions.R")
source("./code/registry_config.R")

default_target_file_dir <- get_registry_file_dir("bt_simu_gp_uq")
if (!exists("target_file_dir", inherits = TRUE)) {
  target_file_dir <- default_target_file_dir
}

if (sys.nframe() == 0 && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    target_file_dir <- args[[1]]
  }
}

target_registry_location <- get_registry_location("bt_simu_gp_uq")
target_registry_location$file.dir <- target_file_dir
source_registry_location <- get_registry_location("bt_simulation")

get_gp_uq_vs_wvl <- function(obj, job) {
  wvl <- as.numeric(rownames(xt))
  x_imp <- obj$x_pred

  x_q <- obj$x_sd * qnorm(0.975)
  S.test.w <- which(apply(obj$x, 2, function(v) all(is.na(v))))
  S.test.o <- which(is.na(obj$x), arr.ind = TRUE)
  S.test.o <- S.test.o[-which(S.test.o[, "col"] %in% S.test.w), ]

  Ome.test <- is.na(obj$x)
  Ome.test.o <- matrix(FALSE, nrow(obj$x), ncol(obj$x))
  Ome.test.o[S.test.o] <- TRUE

  Ome.test.w <- matrix(FALSE, nrow(obj$x), ncol(obj$x))
  Ome.test.w[, S.test.w] <- TRUE
  Ome.test.w[obj$cp$S.test.o2w] <- TRUE

  df <- data.frame(
    wvl = wvl,
    overall = apply((abs(x_imp - xt) <= x_q) * Ome.test, 1, function(v) sum(v)) / apply(Ome.test, 1, sum),
    o = apply((abs(x_imp - xt) <= x_q) * Ome.test.o, 1, function(v) sum(v)) / apply(Ome.test.o, 1, sum),
    w = apply((abs(x_imp - xt) <= x_q) * Ome.test.w, 1, function(v) sum(v)) / apply(Ome.test.w, 1, sum)
  )

  x_q[!is.na(obj$x)] <- NA
  mean_length.w <- apply(x_q[, S.test.w, drop = FALSE], 1, mean)
  x_q[, S.test.w] <- NA
  mean_length.o <- apply(x_q, 1, mean, na.rm = TRUE)

  list(df = df, q.w = mean_length.w, q.o = mean_length.o)
}

extract_job_par_column <- function(pars, name, default = NULL, transform = identity, prototype = default) {
  if (!name %in% names(pars)) {
    return(rep(default, nrow(pars)))
  }

  column <- pars[[name]]
  vapply(seq_len(nrow(pars)), function(i) {
    value <- if (is.list(column)) column[[i]] else column[i]
    transform(value)
  }, FUN.VALUE = prototype)
}

get_target_job_pars <- function(target) {
  pars <- try(unwrap(getJobPars(reg = target)), silent = TRUE)
  if (inherits(pars, "try-error")) { # unwrap failed, likely because there are duplicated column names
    return(data.frame())
  }

  data.frame(
    job.id = extract_job_par_column(pars, "job.id", NA_integer_, as.integer, NA_integer_),
    .id = extract_job_par_column(pars, ".id", NA_integer_, as.integer, NA_integer_)
  )
}

sync_target_job_tags <- function(target, reg) {
  target_pars <- get_target_job_pars(target)
  if (nrow(target_pars) == 0 || !".id" %in% names(target_pars)) {
    return(invisible(NULL))
  }

  raw_parent_pars <- getJobPars(reg = reg)
  parent_pars <- data.frame(
    parent_job.id = extract_job_par_column(raw_parent_pars, "job.id", NA_integer_, as.integer, NA_integer_),
    algorithm = extract_job_par_column(raw_parent_pars, "algorithm", NA_character_, as.character, NA_character_),
    p.dt = extract_job_par_column(raw_parent_pars, "prob.pars", list(p.dt = NA_real_), function(x) as.numeric(x$p.dt), NA_real_),
    parent_diff = extract_job_par_column(raw_parent_pars, "algo.pars", list(diff = FALSE), function(x) isTRUE(x$diff), FALSE)
  ) %>%
    filter(parent_job.id %in% unique(target_pars$.id))

  id_maps <- target_pars %>%
    inner_join(parent_pars, by = c(".id" = "parent_job.id"))

  if (nrow(id_maps) == 0) {
    return(invisible(NULL))
  }

  non_diff_ids <- unique(id_maps$job.id[!id_maps$parent_diff])
  if (length(non_diff_ids) > 0) {
    removeJobTags(ids = non_diff_ids, tags = "diff", reg = target)
  }

  for (i in seq_len(nrow(id_maps))) {
    tags <- c(as.character(id_maps$algorithm[i]), as.character(id_maps$p.dt[i]))
    if (id_maps$parent_diff[i]) {
      tags <- c(tags, "diff")
    }
    addJobTags(ids = id_maps$job.id[i], tags = tags, reg = target)
  }

  invisible(NULL)
}

target <- try(
  loadRegistry(
    work.dir = target_registry_location$work.dir,
    file.dir = target_registry_location$file.dir,
    writeable = TRUE,
    make.default = FALSE
  ),
  silent = TRUE
)

if (inherits(target, "try-error")) {
  target <- makeRegistry(
    work.dir = target_registry_location$work.dir,
    file.dir = target_registry_location$file.dir,
    make.default = FALSE
  )
} else {
  if (!interactive()) {
    stop(sprintf(
      "Target registry already exists at '%s'. Re-run interactively to confirm clearing it.",
      target_file_dir
    ))
  }

  clear_target <- readline(sprintf(
    "Target registry already exists at '%s'. Clear it? [y/N]: ",
    target_file_dir
  ))
  if (!tolower(trimws(clear_target)) %in% c("y", "yes")) {
    stop("Aborted without clearing the existing target registry.")
  }
  clearRegistry(reg = target)
}

dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt
rm(dl)
batchExport(export = list(xt = xt), reg = target)

reg <- loadRegistry(
  work.dir = source_registry_location$work.dir,
  file.dir = source_registry_location$file.dir,
  writeable = FALSE
)

source_ids <- findExperiments(prob.name = "ssi", algo.name = "gp", reg = reg)$job.id
if (length(source_ids) > 0) {
  batchMapResults(get_gp_uq_vs_wvl, ids = source_ids, target = target, source = reg)
}

sync_target_job_tags(target = target, reg = reg)

resources <- list(account = "stats_dept1", walltime = "1:00:00", memory = "1000m", ncpus = 4)
resources$chunks.as.arrayjobs <- TRUE
jobs_per_chunk <- 30
submitted_jobs <- findNotDone(reg = target)$job.id

if (length(submitted_jobs) > 0) {
  njobs <- length(submitted_jobs)
  jobdf <- data.frame(
    job.id = submitted_jobs,
    chunk = rep(seq_len(ceiling(njobs / jobs_per_chunk)), each = jobs_per_chunk)[seq_len(njobs)]
  )
  submitJobs(jobdf, resources = resources, reg = target)
}

cat("GP UQ jobs submitted.\n")
