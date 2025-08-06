# Check if running on local machine or slurm clusters
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))

if (!is.na(cores)) {
  message("Configuring batchtools for SLURM cluster")
  cluster.functions = makeClusterFunctionsSlurm() # configure for SLURM cluster
} else {
  message("Configuring batchtools for LOCAL execution")
  cores <- max(parallel::detectCores() - 1, 1)
  cluster.functions = makeClusterFunctionsMulticore(ncpus = cores) # configure for local machine
}
