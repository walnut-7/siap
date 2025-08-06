# Select maximum time lag for MARSS.
args <- commandArgs(trailingOnly = TRUE)
account = args[1] # slurm account name

# ---------- prepare ---------
library(batchtools)
source("code/read_data_functions.R")
source("code/preprocess_functions.R")
source("code/postprocess_functions.R")
library(dplyr)

# --------- data ----------
dl <- load_ssi(option = "true")
x <- dl$x
wvl <- as.numeric(rownames(x))
t <- ymd(colnames(x))
d1 <- nrow(x)
d2 <- ncol(x)

pre_par <- preprocess_func(x)
x_ready <- pre_par$x_ready

# ----------- ward clustering -----------
x_cmplt <- x_ready[, -which(apply(x_ready, 2, function(v) any(is.na(v))))] # x_ready with complete columns

dist_matrix <- dist(x_cmplt) # Euclidian distance

hc <- hclust(dist_matrix, method = "ward.D2") # Perform Ward clustering

# plot(hc)

ngp = 10
group <- cutree(hc, k = ngp)

# ------------- acf --------------
group_centers <- sapply(1:ngp, function(i) apply(x_cmplt[which(group == i),], 2, mean)) %>% t() # group centers
acf <- apply(group_centers, 1, function(v) acf(v, lag.max = 100)$acf) %>% t()

# ------------- select time lag ------------------
lags <- apply(acf, 1, function(v) max(which(v >= 0.2))) # auto-correlation threshold: 0.2
saveRDS(max(lags), file = "./output/realdata/marss_p.rds")
