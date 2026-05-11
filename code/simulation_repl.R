args <- commandArgs(trailingOnly = TRUE)
# cfg_path <- args[1]
pdt <- as.numeric(args[1])
repl <- as.integer(args[2])
use_calibration <- if (length(args) >= 3) {
  as.logical(as.integer(args[3]))
} else {
  FALSE
}
out_path <- paste0("./code/repl/mask_pdt", pdt*10, '_repl', repl, ".json")
library(jsonlite)
library(batchtools)
source("code/read_data_functions.R")

seeds <- repl
d1 <- 2104
d2 <- 1783

# subsample = function(d1, d2, p.dt, p.sctt, cal = F, ...) {
#   mask <- chunk_mask(m = d1, n = d2, option = "uniform", ratio = p.dt, return = "mask")
#   ind <- which(rbinom(d1*d2, 1, p.sctt) == 1)
#   mask[ind] <- T
#   return(which(mask == T)-1) # return missing entries column-major indices, starting from 0.
# }

subsample_cal = function(d1, d2, p.dt, p.sctt, cal=F, p.cal = 0.1, ...) {
  mask <- chunk_mask(m = d1, n = d2, option = "uniform", ratio = p.dt, return = "mask")
  ind <- which(rbinom(d1*d2, 1, p.sctt) == 1)
  mask[ind] <- T
  
  # ---------- cp -----------
  if (cal == T) {
    S_full <- which(apply(mask, 2, function(v) all(v == T))) # fully missing
    S.cal.w <- setdiff(1:ncol(mask), S_full)[as.logical(rbinom(ncol(mask) - length(S_full), 1, p.cal))]
    M.o <- (mask==F)
    rownames(M.o) <- 1:nrow(mask)
    M.o[, S.cal.w] <- F
    S.cal.o <- which(M.o == T, arr.ind = T)[as.logical(rbinom(sum(M.o), 1, p.cal)), ]
    return(list(miss = which(mask == T)-1, cp = list(S.cal.o = S.cal.o-1, S.cal.w = S.cal.w-1))) 
  }
  return(which(mask == T)-1) # return missing entries
}


missing <- vector('list', length(seeds))
for (i in 1:length(seeds)) {
  seed <- seeds[i]
  set.seed(seed)
  mask <- subsample_cal(d1, d2, p.dt = pdt, p.sctt = pdt, cal = use_calibration)
  missing[[i]] <- mask
}
write_json(missing, out_path, digits = NA)
print(paste0('Repl ', repl, ' with pdt ', pdt, ' generated (calibration ', use_calibration, ').'))
