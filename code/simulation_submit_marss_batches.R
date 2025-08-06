args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

account = ... # change to your slurm account name

library(batchtools)
reg = loadRegistry(work.dir = getwd(), file.dir = "bt_simulation", writeable = T)
marss_ids <- findExperiments(prob.name = "ssi_batch", algo.pars = (model == "low-rank"))$job.id
submit_ids = marss_ids

resources = list(account = account, walltime = '10:00:00', memory='1000m', ncpus=4) # MARSS
resources$chunks.as.arrayjobs = TRUE 

jobs_per_chunk = 10 
njobs <- length(submit_ids)
jobdf <- data.frame(job.id = submit_ids, 
                    chunk=rep(1:ceiling(njobs / jobs_per_chunk), jobs_per_chunk)[1:njobs],
                    submit_batch=rep(1:ceiling(njobs / 4000), each = 4000)[1:njobs])

submitJobs(jobdf[which(jobdf$submit_batch == i), ], resources = resources)
