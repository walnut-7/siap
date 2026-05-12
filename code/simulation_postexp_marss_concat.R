# Usage:
#   Rscript code/simulation_postexp_marss_concat.R [target_file_dir]
#   target_file_dir <- "/path/to/bt_simu_marss"; source("code/simulation_postexp_marss_concat.R")

library(dplyr)
library(batchtools)

source("./code/read_data_functions.R")
source("./code/registry_config.R")

default_target_file_dir <- get_registry_file_dir("bt_simu_marss")
if (!exists("target_file_dir", inherits = TRUE)) {
  target_file_dir <- default_target_file_dir
}

if (sys.nframe() == 0 && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    target_file_dir <- args[[1]]
  }
}

target_registry_location <- get_registry_location("bt_simu_marss")
target_registry_location$file.dir <- target_file_dir
source_registry_location <- get_registry_location("bt_simulation")

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

extract_job_par_list_column <- function(pars, name, default = integer()) {
  if (!name %in% names(pars)) {
    return(rep(list(default), nrow(pars)))
  }

  column <- pars[[name]]
  lapply(seq_len(nrow(pars)), function(i) {
    value <- if (is.list(column)) column[[i]] else column[i]
    if (is.null(value)) default else value
  })
}

build_source_job_map <- function(reg) {
  source_ids <- findExperiments(
    prob.name = "ssi_batch",
    algo.name = "marss",
    algo.pars = (model == "low-rank"),
    reg = reg
  )$job.id

  if (length(source_ids) == 0) {
    return(data.frame())
  }

  raw_pars <- getJobTable(source_ids, reg = reg)
  data.frame(
    source_job.id = extract_job_par_column(raw_pars, "job.id", NA_integer_, as.integer, NA_integer_),
    repl = extract_job_par_column(raw_pars, "repl", NA_integer_, as.integer, NA_integer_),
    pdt = extract_job_par_column(
      raw_pars,
      "prob.pars",
      list(p.dt = NA_real_),
      function(x) as.numeric(x$p.dt),
      NA_real_
    ),
    block_id = extract_job_par_column(
      raw_pars,
      "prob.pars",
      list(id = NA_integer_),
      function(x) as.integer(x$id),
      NA_integer_
    ),
    parent_diff = extract_job_par_column(
      raw_pars,
      "algo.pars",
      list(diff = FALSE),
      function(x) isTRUE(x$diff),
      FALSE
    )
  ) %>%
    arrange(repl, pdt, parent_diff, block_id, source_job.id)
}

build_concat_args <- function(source_job_map) {
  if (nrow(source_job_map) == 0) {
    return(data.frame())
  }

  grouped <- split(
    source_job_map,
    list(source_job_map$repl, source_job_map$pdt, source_job_map$parent_diff),
    drop = TRUE
  )

  do.call(
    rbind,
    lapply(grouped, function(df) {
      df <- df[order(df$block_id, df$source_job.id), , drop = FALSE]
      data.frame(
        repl = df$repl[1],
        pdt = df$pdt[1],
        diff_flag = df$parent_diff[1],
        source_ids = I(list(df$source_job.id)),
        stringsAsFactors = FALSE
      )
    })
  ) %>%
    arrange(repl, pdt, diff_flag)
}

get_target_job_layout <- function(target) {
  pars <- try(getJobPars(reg = target), silent = TRUE)
  if (inherits(pars, "try-error") || nrow(pars) == 0) {
    return(data.frame())
  }

  layout <- data.frame(
    repl = extract_job_par_column(pars, "repl", NA_integer_, as.integer, NA_integer_),
    pdt = extract_job_par_column(pars, "pdt", NA_real_, as.numeric, NA_real_),
    diff_flag = extract_job_par_column(pars, "diff_flag", NA, isTRUE, NA),
    stringsAsFactors = FALSE
  )
  layout$source_ids_key <- vapply(
    extract_job_par_list_column(pars, "source_ids"),
    function(ids) paste(as.integer(ids), collapse = ","),
    FUN.VALUE = character(1)
  )

  layout %>% arrange(repl, pdt, diff_flag, source_ids_key)
}

same_target_layout <- function(target, concat_args) {
  target_layout <- get_target_job_layout(target)
  if (nrow(target_layout) != nrow(concat_args)) {
    return(FALSE)
  }

  expected_layout <- concat_args %>%
    mutate(source_ids_key = vapply(source_ids, function(ids) paste(as.integer(ids), collapse = ","), FUN.VALUE = character(1))) %>%
    select(repl, pdt, diff_flag, source_ids_key) %>%
    arrange(repl, pdt, diff_flag, source_ids_key)

  identical(target_layout, expected_layout)
}

sync_target_job_tags <- function(target, reg) {
  target_pars <- try(unwrap(getJobPars(reg = target)), silent = TRUE)
  if (inherits(target_pars, "try-error") || nrow(target_pars) == 0) {
    return(invisible(NULL))
  }

  target_map <- data.frame(
    job.id = extract_job_par_column(target_pars, "job.id", NA_integer_, as.integer, NA_integer_),
    pdt = extract_job_par_column(target_pars, "pdt", NA_real_, as.numeric, NA_real_)
  )
  target_map$source_ids <- extract_job_par_list_column(target_pars, "source_ids")

  parent_ids <- unique(unlist(target_map$source_ids, use.names = FALSE))
  raw_parent_pars <- getJobPars(parent_ids, reg = reg)
  parent_pars <- data.frame(
    parent_job.id = extract_job_par_column(raw_parent_pars, "job.id", NA_integer_, as.integer, NA_integer_),
    parent_diff = extract_job_par_column(
      raw_parent_pars,
      "algo.pars",
      list(diff = FALSE),
      function(x) isTRUE(x$diff),
      FALSE
    )
  )

  parent_diff_lookup <- setNames(parent_pars$parent_diff, parent_pars$parent_job.id)
  target_diff <- vapply(target_map$source_ids, function(ids) {
    parent_diff <- unique(parent_diff_lookup[as.character(ids)])
    parent_diff <- parent_diff[!is.na(parent_diff)]
    if (length(parent_diff) == 0) {
      FALSE
    } else {
      if (length(parent_diff) > 1) {
        stop("Target job mixes diff and non-diff MARSS parents.")
      }
      parent_diff[[1]]
    }
  }, FUN.VALUE = FALSE)

  non_diff_ids <- target_map$job.id[!target_diff]
  if (length(non_diff_ids) > 0) {
    removeJobTags(ids = non_diff_ids, tags = "diff", reg = target)
  }

  diff_ids <- target_map$job.id[target_diff]
  if (length(diff_ids) > 0) {
    addJobTags(ids = diff_ids, tags = "diff", reg = target)
  }

  invisible(NULL)
}

concat_marss <- function(repl, pdt, diff_flag, source_ids) {
  library(batchtools)

  reg <- loadRegistry(
    work.dir = source_registry_location$work.dir,
    file.dir = source_registry_location$file.dir,
    writeable = FALSE
  )

  fn <- function(obj1, obj2) {
    list(
      x = rbind(obj1$x, obj2$x),
      x_pred = rbind(obj1$x_pred, obj2$x_pred),
      x_sd = rbind(obj1$x_sd, obj2$x_sd),
      rae = list(
        rae.w = rbind(obj1$rae$rae.w, obj2$rae$rae.w),
        rae.o = rbind(obj1$rae$rae.o, obj2$rae$rae.o)
      )
    )
  }

  reduceResults(fn, ids = source_ids, reg = reg)
}

dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt
rm(dl)

reg <- loadRegistry(
  work.dir = source_registry_location$work.dir,
  file.dir = source_registry_location$file.dir,
  writeable = FALSE
)

source_job_map <- build_source_job_map(reg)
if (nrow(source_job_map) == 0) {
  stop("No MARSS source jobs found in the source registry.")
}

concat_args <- build_concat_args(source_job_map)

target <- try(
  loadRegistry(
    work.dir = target_registry_location$work.dir,
    file.dir = target_registry_location$file.dir,
    writeable = TRUE,
    make.default = FALSE
  ),
  silent = TRUE
)

remap_target <- FALSE
if (inherits(target, "try-error")) {
  target <- makeRegistry(
    work.dir = target_registry_location$work.dir,
    file.dir = target_registry_location$file.dir,
    seed = 1,
    make.default = FALSE
  )
  remap_target <- TRUE
} else if (!same_target_layout(target, concat_args)) {
  if (!interactive()) {
    stop(sprintf(
      "Target registry layout has changed for '%s'. Re-run interactively to confirm clearing the existing target registry.",
      target_file_dir
    ))
  }

  clear_target <- readline(sprintf(
    "Target registry layout has changed. Clear the existing target registry at '%s'? [y/N]: ",
    target_file_dir
  ))
  if (!tolower(trimws(clear_target)) %in% c("y", "yes")) {
    stop("Aborted without clearing the existing bt_simu_marss registry.")
  }

  clearRegistry(reg = target)
  remap_target <- TRUE
}

if (remap_target) {
  batchExport(export = list(xt = xt), reg = target)
  batchMap(fun = concat_marss, args = concat_args, reg = target)
}

sync_target_job_tags(target = target, reg = reg)

resources <- list(account = "stats_dept1", walltime = "1:00:00", memory = "1000m", ncpus = 4)
resources$chunks.as.arrayjobs <- TRUE
jobs_per_chunk <- 20
submitted_jobs <- findNotDone(reg = target)$job.id

if (length(submitted_jobs) > 0) {
  njobs <- length(submitted_jobs)
  jobdf <- data.frame(
    job.id = submitted_jobs,
    chunk = rep(seq_len(ceiling(njobs / jobs_per_chunk)), each = jobs_per_chunk)[seq_len(njobs)]
  )
  submitJobs(reg = target, resources = resources, ids = jobdf)
}

cat("MARSS concatenation jobs submitted.\n")
