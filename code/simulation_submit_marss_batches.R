args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
use_diff <- ifelse(length(args) >= 2, as.logical(args[2]), TRUE)
snapshot_path <- ifelse(length(args) >= 3, args[3], NA_character_)
diff_mode <- ifelse(isTRUE(use_diff), "diff-only", "all-pending")
cat(sprintf("[submit-marss] init | submit_batch=%s | mode=%s | snapshot=%s\n",
            i, diff_mode, snapshot_path))

library(batchtools)
library(dplyr)
reg = loadRegistry(work.dir = "/home/yuxuank/solar", file.dir = "/nfs/turbo/lsa-ychenang/yuxuank/bt_simulation", writeable = T)

if (!is.na(snapshot_path) && file.exists(snapshot_path)) {
  marss_ids <- scan(snapshot_path, what = integer(), quiet = TRUE)
  cat(sprintf("[submit-marss] loaded %d job ids from snapshot.\n", length(marss_ids)))
} else {
  if (isTRUE(use_diff)) {
    marss_ids <- (
      findExperiments(prob.name = "ssi_batch",algo.pars = (model == "low-rank")) %>%
        ijoin(findTagged("diff")) %>%
        ijoin(findNotStarted())
    )$job.id
  } else {
    marss_ids <- (
      findExperiments(prob.name = "ssi_batch", algo.pars = (model == "low-rank" )) %>%
        ijoin(findNotStarted())
    )$job.id
  }

  if (!is.na(snapshot_path)) {
    writeLines(as.character(marss_ids), con = snapshot_path)
    cat(sprintf("[submit-marss] created snapshot with %d job ids.\n", length(marss_ids)))
  }
}
submit_ids = marss_ids

resources = list(account = "xianglei1", walltime = '10:00:00', memory='1000m', ncpus=4) # MARSS
resources$chunks.as.arrayjobs = TRUE 

jobs_per_chunk = 10 
njobs <- length(submit_ids)
cat(sprintf("[submit-marss] discovery | pending_jobs=%d | mode=%s\n", njobs, diff_mode))
if (njobs == 0) {
  cat("[submit-marss] no jobs to submit; exiting.\n")
  quit(save = "no", status = 0)
}

submit_batches <- ceiling(njobs / 4000)

jobdf <- data.frame(job.id = submit_ids, 
                    chunk=rep(1:ceiling(njobs / jobs_per_chunk), jobs_per_chunk)[1:njobs],
                    submit_batch=rep(1:ceiling(njobs / 4000), each = 4000)[1:njobs])
selected_jobs <- jobdf[which(jobdf$submit_batch == i), , drop = FALSE]
if (nrow(selected_jobs) == 0) {
  cat(sprintf("[submit-marss] batch %d exceeds available batches (%d); exiting.\n",
              i, submit_batches))
  quit(save = "no", status = 0)
}
cat(sprintf("[submit-marss] submitting batch %d/%d: %d jobs across %d chunks\n",
            i, submit_batches, nrow(selected_jobs), length(unique(selected_jobs$chunk))))
cat(sprintf("[submit-marss] job.id range: %d - %d\n",
            min(selected_jobs$job.id), max(selected_jobs$job.id)))
if (nrow(selected_jobs) <= 20) {
  cat(sprintf("[submit-marss] job.id list: %s\n", paste(selected_jobs$job.id, collapse = ", ")))
}

cat("[submit-marss] dispatching jobs to scheduler...\n")
submitJobs(selected_jobs, resources = resources)
cat("[submit-marss] submitJobs completed successfully.\n")
