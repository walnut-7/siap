source("code/read_data_functions.R")
source("code/plot_functions.R")

library(dplyr)
library(ggplot2)
library(ggbreak)
library(scales)
library(ggpubr)
library(lubridate)
library(tidyverse)

dl <- load_ssi(option = "true")
x <- dl$x
wvl <- as.numeric(rownames(x))
t <- ymd(colnames(x))
d1 <- nrow(x)
d2 <- ncol(x)

# ------------ Figures -------------
## ------- load results ----------
siap.res<- readRDS(file = "./output/realdata/siap_traintest.rds")
gp.res<- readRDS(file = "./output/realdata/gp_traintest.rds")
marss.res<- readRDS(file = "./output/realdata/marss_lowrank_traintest.rds")

### ------ siap --------
siap_imp <- siap.res$siap$fit2$x_imp
colnames(siap_imp) <- colnames(x)
rownames(siap_imp) <- wvl

siap_q <- matrix(0, nrow = nrow(x), ncol = ncol(x))
for (i in 1:nrow(siap.res$cp$S.test.o2w)) {
  siap_q[siap.res$cp$S.test.o2w[i,,drop=F]] <- siap.res$cp$cp_q.w[siap.res$cp$S.test.o2w[i,1]]
}
S.test.w <- which(apply(siap.res$x, 2, function(v) all(is.na(v))))
siap_q[,S.test.w] <- siap.res$cp$cp_q.w

S.cal.w <- siap.res$cp$S.cal.w
S.test.o <- which(is.na(siap.res$x), arr.ind = T)
S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.test.w), ]
S.test.o <- S.test.o[-which(S.test.o[,"col"] %in% S.cal.w), ]
for (i in 1:nrow(S.test.o)) {
  siap_q[S.test.o[i,,drop=F]] <- siap.res$cp$cp_q.o[S.test.o[i,1]]
}
colnames(siap_q) <- colnames(x)
rownames(siap_q) <- wvl

### ------- gp --------
gp_imp <- gp.res$x_pred
gp_sims <- gp.res$x_sims
gp_sd <- gp.res$x_sd

colnames(gp_imp) <- colnames(x)
rownames(gp_imp) <- wvl
colnames(gp_sims) <- colnames(x)
rownames(gp_sims) <- wvl

### ------ marss -------
marss_imp <-marss.res$x_pred
marss_boots <- marss.res$x_boots
marss_sd <- marss.res$x_sd

colnames(marss_imp) <- colnames(x)
rownames(marss_imp) <- wvl
colnames(marss_boots) <- colnames(x)
rownames(marss_boots) <- wvl

### --------- softimpute model ---------
siap.res<- readRDS(file = "./output/realdata/siap_traintest.rds")
x_c <- (siap.res$x - apply(siap.res$x, 1, mean, na.rm = T)) / apply(siap.res$x, 1, sd, na.rm = T)
start <- Sys.time()
si.res <- softImpute::softImpute(x_c, rank.max = 10, lambda = 1)
si_runtime <- Sys.time() - start
si_imp <- (si.res$u %*% (t(si.res$v) * si.res$d)) * apply(siap.res$x, 1, sd, na.rm = T) + apply(siap.res$x, 1, mean, na.rm = T)

rr <- abs((si_imp - x)/x)
si_rae.w <- rr[,siap.res$S.test.w,drop=F]

mask.o <- matrix(NA, d1, d2)
mask.o[siap.res$S.test.o] <- 1
si_rae.o <- rr*mask.o

colnames(si_imp) <- colnames(x)
rownames(si_imp) <- wvl

## ----------- rae margin to naive mean -------------
{
  baseline = matrix(apply(x, 1, mean, na.rm=T), nrow(x), ncol(x))
  baseline = abs((baseline - x)/baseline)
  bl.w <- baseline[,siap.res$S.test.w,drop=F]
  mask.o <- matrix(NA, nrow(x), ncol(x))
  mask.o[siap.res$S.test.o] <- 1
  bl.o <- baseline*mask.o
  
  df_rae <- data.frame(siap = c(bl.w)-c(siap.res$rae$rae.w),
                       gp = c(bl.w)-c(gp.res$rae$rae.w),
                       marss = c(bl.w)-c(marss.res$rae$rae.w),
                       si = c(bl.w) - c(si_rae.w),
                       type = "w")
  df_rae <- rbind(df_rae, data.frame(siap = c(bl.o)-c(siap.res$rae$rae.o),
                                     gp = c(bl.o)-c(gp.res$rae$rae.o),
                                     marss = c(bl.o)-c(marss.res$rae$rae.o),
                                     si = c(bl.o) - c(si_rae.o),
                                     type = "o"))
  rae_stats <- df_rae %>% 
    group_by(type) %>%
    summarise(
      siap_mean = mean(siap, na.rm = TRUE),
      siap_q = sd(siap, na.rm = TRUE)/sqrt(sum(!is.na(siap))),
      gp_mean = mean(gp, na.rm = TRUE),
      gp_q = sd(gp, na.rm = TRUE)/sqrt(sum(!is.na(gp))),
      marss_mean = mean(marss, na.rm = TRUE),
      marss_q = sd(marss, na.rm = TRUE)/sqrt(sum(!is.na(marss))),
      si_mean = mean(si, na.rm = TRUE),
      si_q = sd(si, na.rm = TRUE)/sqrt(sum(!is.na(si)))
    )
  
  long_stats <- rae_stats %>%
    pivot_longer(
      cols = starts_with("siap") | starts_with("gp") | starts_with("marss") | starts_with("si"),
      names_to = c("method", ".value"),
      names_pattern = "(.*)_(.*)"
    )
  
  p <- ggplot(long_stats, aes(x = method, y = mean, col = type)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.5) +  # Scatter points
    geom_errorbar(aes(ymin = mean - q, ymax = mean + q),
                  width = 0.2, position = position_dodge(width = 0.5)) +  # Error bars
    labs(title = "", x = "Method", y = "MRAE margin w.r.t. naive mean") +
    scale_color_manual(name = "Type", labels = c("scattered", "downtime"), values = c(2, 3)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")    # Place the legend at the top
  
  ggsave(paste0("./output/realdata/mrae_margin_to_mean_imputation.jpeg"), plot = p, width = 6, height = 4, dpi = 300)
}

# ------------ Table -------------
## ----------- runtime -------------
{
  siap.res$siap$fit1$time + siap.res$siap$fit2$time
  gp.res$runtime$runtime.gp + gp.res$runtime$runtime.pred + gp.res$runtime$runtime.sims
  as.numeric(sum(marss.res$runtime$runtime.marss) + sum(marss.res$runtime$runtime.boot), units = "mins")
  as.numeric(si_runtime, units = "mins")
}

## ---------- coverage --------------
{
  coverage_df <- data.frame(siap = c(siap.res$ave_coverage, siap.res$ave_coverage.w, siap.res$ave_coverage.o),
                            gp = c(gp.res$ave_coverage, gp.res$ave_coverage.w, gp.res$ave_coverage.o),
                            marss = c((marss.res$ttl_coverage.w + marss.res$ttl_coverage.o) / (marss.res$n.w + marss.res$n.o),
                                      marss.res$ttl_coverage.w/marss.res$n.w,
                                      marss.res$ttl_coverage.o/marss.res$n.o))
  
  rownames(coverage_df) <- c("overall", "downtime", "scattered")
  write.table(coverage_df, file = "./output/realdata/coverage_traintest.csv", sep = "\t")
}

## ------------- time series of a row ------------
{
  i = 200
  x_coordinates = which(apply(siap.res$x, 2, function(t) all(is.na(t))))
  siap_imp1 = siap_imp
  siap_imp1[!is.na(siap.res$x)] = x[!is.na(siap.res$x)] # fill calibration set with observed values
  jpeg(file=paste0("./output/realdata/timeseries.jpeg"), width=10, height=3, units = 'in', res = 300, pointsize=12)
  plot(t, x[i,], type = "n", lty = 2, ylab = paste0("SSI, ", i, "-th row (", round(wvl[i],2), "nm)"), xlab = "Day")
  for (j in seq_along(x_coordinates)) {
    rect(xleft = t[x_coordinates[j]] - 0.5, xright = t[x_coordinates[j]] + 0.5, ybottom = min(x[i,],na.rm = T), ytop = max(x[i,],na.rm = T), col = "gray", border = NA)
  } # downtime
  lines(t,si_imp[i,], col = "violet")
  lines(t,siap_imp1[i,], col = "red")
  lines(t,marss_imp[i,], col = "lightgreen") 
  lines(t,gp_imp[i,], col = "lightblue")
  lines(t,x[i,])
  legend("topleft", lty = 1, col = c("red", "lightgreen", "lightblue", "violet"), legend = c("SIAP", "MARSS", "GP", "SoftImpute-ALS"))
  dev.off()
}

## ---------- snapshot of reconstruction result -------------
{  
  zlim = range(c(x,
                 siap.res$siap$fit2$x_imp,
                 gp.res$x_pred,
                 marss.res$x_pred), na.rm=T)
  
  jpeg(file=paste0("./output/realdata/snapshot.jpeg"), width=20, height=12, units = 'in', res = 100, pointsize=16)
  par(mfrow = c(2, 2), mar = c(3,5,5,8))
  panel_image(x, yleg = wvl, xlab = "", xleg = t, title = "TSIS-1 SSI Observations", zlim = zlim, legend.strip = F)
  panel_image(siap.res$siap$fit2$x_imp, yleg = wvl, xlab = "", xleg = t, title = "SIAP", zlim = zlim, legend.strip = F)
  panel_image(gp.res$x_pred, yleg = wvl, xlab = "", xleg = t, title = "GP", zlim = zlim, legend.strip = F)
  panel_image(marss.res$x_pred, yleg = wvl, xlab = "", xleg = t, title = "MARSS", zlim = zlim, legend.strip = F)
  image.plot(legend.only = TRUE, zlim = zlim, col = tim.colors(), legend.args = list(text = expression("Irradiance" ~ "[" * mW ~ m^{-2} * nm^{-1} * "]" ), side = 1, line = 2), legend.mar = 1, horizontal = F)
  #dev.copy(png, file=paste0("./test9/example_", 10*pdt, ".png"), width=20, height=20, units = 'in', res = 300, pointsize=10)
  dev.off()
}

## --------- snapshot of half interval width ---------------
{
  gp_q <- gp.res$x_sd * qnorm(0.975)
  marss_q <- marss.res$x_sd * qnorm(0.975)
  
  jpeg(file=paste0("./output/realdata/snapshot_q.jpeg"), width=15, height=4, units = 'in', res = 300, pointsize=12)
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
  par(mar = c(3,4,4,8))
  panel_image(siap_q, yleg = wvl, xlab = "", xleg = t, title = "SIAP", legend = expression("[" * mW ~ m^{-2} * nm^{-1} * "]" ), legend.mar = 5, legend.side = 1)
  panel_image(gp_q, yleg = wvl, xlab = "", xleg = t, title = "GP", legend = expression("[" * mW ~ m^{-2} * nm^{-1} * "]" ), legend.mar = 5, legend.side = 1)
  #par(mar = c(3,4,4,8))
  panel_image(marss_q, yleg = wvl, xlab = "", xleg = t, title = "MARSS", legend = expression("[" * mW ~ m^{-2} * nm^{-1} * "]" ), legend.mar = 5, legend.side = 1)
  #image.plot(legend.only = TRUE, zlim = range(marss_q), col = tim.colors(), legend.args = list(text = expression("[" * mW ~ m^{-2} * nm^{-1} * "]" ), side = 1, line = 2), legend.mar = 1, horizontal = F)
  dev.off()
}

## --------- snapshot of sd ---------------
{
  gp_sd <- gp.res$x_sd
  marss_sd <- marss.res$x_sd
  
  zlim = range(c(siap_q/qnorm(0.975), gp_sd, marss_sd))
  
  jpeg(file=paste0("./output/realdata/snapshot_sd.jpeg"), width=15, height=4, units = 'in', res = 300, pointsize=12)
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1.3))
  par(mar = c(3,4,4,1))
  panel_image(siap_q/qnorm(0.975), yleg = wvl, xlab = "", xleg = t, title = "SIAP", zlim = zlim, legend.strip = F)
  panel_image(gp_sd, yleg = wvl, xlab = "", xleg = t, title = "GP", zlim = zlim, legend.strip = F)
  par(mar = c(3,4,4,8))
  panel_image(marss_sd, yleg = wvl, xlab = "", xleg = t, title = "MARSS", zlim = zlim, legend.strip = F)
  image.plot(legend.only = TRUE, zlim = zlim, col = tim.colors(), legend.args = list(text = expression("[" * mW ~ m^{-2} * nm^{-1} * "]" ), side = 1, line = 2), legend.mar = 1, horizontal = F)
  dev.off()
}

## ------ daily relative percentage difference in test set -------
### downtime
{
  j = 724 # 2020-03-06
  
  r.siap = (siap_imp[,j] - x[,j])/x[,j]*100
  r.gp = (gp_imp[,j] - x[,j])/x[,j]*100
  r.marss = (marss_imp[,j] - x[,j])/x[,j]*100
  
  df_low = data.frame(
    wvl = wvl[wvl<400],
    siap = r.siap[wvl<400],
    gp = r.gp[wvl<400],
    marss = r.marss[wvl<400]
  )
  
  df_low <- pivot_longer(df_low, cols = c(siap, gp, marss), 
                         names_to = "method", 
                         values_to = "r")
  p <- ggplot(df_low, aes(x = wvl, y = r, color = method, fill = method)) +
    geom_line() +  # Line for mean values
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y=0
    labs(x = "Wavelength[nm]", y = "[%]", title = paste0("Relative % Difference from TSIS, ", t[j])) +
    theme_minimal() +
    ylim(-1, 1)
  ggsave(paste0("./output/realdata/test_", t[j], "_low.jpeg"), plot = p, width = 6, height = 2, dpi = 300)
  
  df_high = data.frame(
    wvl = wvl[wvl>=400],
    siap = r.siap[wvl>=400],
    gp = r.gp[wvl>=400],
    marss = r.marss[wvl>=400]
  )
  
  df_high <- pivot_longer(df_high, cols = c(siap, gp, marss), 
                          names_to = "method", 
                          values_to = "r")
  p <- ggplot(df_high, aes(x = wvl, y = r, color = method, fill = method)) +
    geom_line() +  # Line for mean values
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y=0
    labs(x = "Wavelength[nm]", y = "[%]", title ="") +
    theme_minimal() +
    ylim(-0.5, 0.5)
  
  ggsave(paste0("./output/realdata/test_", t[j], "_high.jpeg"), plot = p, width = 6, height = 2, dpi = 300)
}

### scattered
{
  j = siap.res$S.test.o[1,"col"] # 2018-03-14
  
  r.siap = (siap_imp[,j] - x[,j])/x[,j]*100
  r.gp = (gp_imp[,j] - x[,j])/x[,j]*100
  r.marss = (marss_imp[,j] - x[,j])/x[,j]*100
  
  df_low = data.frame(
    wvl = wvl[wvl<400],
    siap = r.siap[wvl<400],
    gp = r.gp[wvl<400],
    marss = r.marss[wvl<400]
  )
  
  df_low <- pivot_longer(df_low, cols = c(siap, gp, marss), 
                         names_to = "method", 
                         values_to = "r")
  p <- ggplot(df_low, aes(x = wvl, y = r, color = method, fill = method)) +
    geom_line() +  # Line for mean values
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y=0
    labs(x = "Wavelength[nm]", y = "[%]", title = paste0("Relative % Difference from TSIS, ", t[j])) +
    theme_minimal() +
    ylim(-0.5, 0.5)
  ggsave(paste0("./output/realdata/test_", t[j], "_low.jpeg"), plot = p, width = 6, height = 2, dpi = 300)
  
  df_high = data.frame(
    wvl = wvl[wvl>=400],
    siap = r.siap[wvl>=400],
    gp = r.gp[wvl>=400],
    marss = r.marss[wvl>=400]
  )
  
  df_high <- pivot_longer(df_high, cols = c(siap, gp, marss), 
                          names_to = "method", 
                          values_to = "r")
  p <- ggplot(df_high, aes(x = wvl, y = r, color = method, fill = method)) +
    geom_line() +  # Line for mean values
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y=0
    labs(x = "Wavelength[nm]", y = "[%]", title = "") +
    theme_minimal() +
    ylim(-0.5, 0.5)
  
  ggsave(paste0("./output/realdata/test_", t[j], "_high.jpeg"), plot = p, width = 6, height = 2, dpi = 300)
}

## --------- compare with csim ssi ------------
### helpers function
{
  integrate_ssi <- function(ssi) {
    # - ssi: Matrix. Row names must be the wavelength
    # No missing values should be in the input data
    wvl <- as.numeric(rownames(ssi))
    area <- colSums((ssi[-1,,drop=F] + ssi[-nrow(ssi),,drop=F]) * (wvl[-1] - wvl[-length(wvl)]) / 2, )
    area <- matrix(area, 1)
    rownames(area) <- paste0(round(range(wvl),1), collapse = "~")
    colnames(area) <- colnames(ssi)
    return(area)
  }
  
  binned <- function(ssi) {
    wvl <- as.numeric(rownames(ssi))
    rbind(
      integrate_ssi(ssi[which(wvl>210)[1]:which(wvl>300)[1],]) / 1000,
      integrate_ssi(ssi[which(wvl>300)[1]:which(wvl>400)[1],]) / 1000,
      integrate_ssi(ssi[which(wvl>700)[1]:which(wvl>1000)[1],]) / 1000,
      integrate_ssi(ssi[which(wvl>1000)[1]:which(wvl>1300)[1],]) / 1000
    )
  }
  
}

{
  csim = load_ssi(option = "csim")
  csim = csim$x * 1000
  t_csim = ymd(colnames(csim))
  t_csim <- as.character(t_csim)
  wvl_csim = as.numeric(rownames(csim))
  
  csim_uniq <- matrix(NA, nrow = nrow(csim), ncol = length(unique(t_csim)))
  colnames(csim_uniq) <- unique(t_csim)
  rownames(csim_uniq) <- wvl_csim
  for (i in unique(t_csim)) {
    csim_uniq[, i] <- rowMeans(csim[, which(t_csim == i), drop = F], na.rm = T)
  }
  
  csim_match <- matrix(NA, nrow = length(wvl), ncol = ncol(csim_uniq))
  for (j in 1:ncol(csim_uniq)) {
    csim_match[,j] <- approx(wvl_csim, csim_uniq[,j], wvl)$y # linear interpolation
  }
  
  t_csim <- unique(t_csim)
  colnames(csim_match) <- as.character(t_csim)
  rownames(csim_match) <- wvl
}

### binned ssi
{
  siap_imp1 = siap_imp
  siap_imp1[!is.na(siap.res$x)] = x[!is.na(siap.res$x)] # fill calibration set with observed values
  binned_siap <- binned(siap_imp1)
  
  binned_siap_q <- sqrt(binned((siap_q/qnorm(0.975))^2)) * qnorm(0.975)
  
  binned_gp <- binned(gp_imp)
  
  x_sims <- gp.res$x_sims
  rownames(x_sims) <- wvl
  binned_gp_sims <- array(dim=c(4, ncol(x), dim(x_sims)[3]))
  for (i in 1:dim(x_sims)[3]) {
    binned_gp_sims[,,i] <- binned(x_sims[,,i])
  }
  binned_gp_intervals <- apply(binned_gp_sims, c(1,2), quantile, probs = c(0.025,0.975), na.rm = T)
  dimnames(binned_gp_intervals) <- list(
    rows = c("low", "high"),
    cols = c("210~300.1", "300.1~400.2",  "701.1~1002.9", "1002.9~1304.2"),
    slices = colnames(x)
  )
  
  gp_q <- gp.res$x_sd * qnorm(0.975)
  marss_q <- marss.res$x_sd * qnorm(0.975)
  
  rownames(gp_q) <- wvl
  binned_gp_q <- binned(gp_q)
  
  binned_marss <- binned(marss_imp)
  
  rownames(marss_q) <- wvl
  binned_marss_q <- binned(marss_q) # conservative
  
  binned_tsis <- binned(x)
  
  binned_csim <- matrix(NA,4,ncol(x))
  rownames(binned_csim) <- c("210~300.1", "300.1~400.2",  "701.1~1002.9", "1002.9~1304.2")
  colnames(binned_csim) <- colnames(x)
  binned_csim[,t_csim] <- binned(csim_uniq)
}

### scale csim
{
  wvl_ranges_ls <- list(c(210,300), c(300,400), c(700,1000), c(1000,1300))
  csim_scale <- rep(NA, 4)
  
  for (i in 1:length(wvl_ranges_ls)) {
    wvl_range <- wvl_ranges_ls[[i]]
    csim_irradiance <- integrate_ssi(csim_uniq[which(wvl_csim>wvl_range[1])[1]:which(wvl_csim>wvl_range[2])[1],]) / 1000
    tsis_irradiance <- binned_tsis[i,t_csim,drop=F]
    j <- which(!is.na(csim_irradiance[40:length(t_csim)])&!is.na(tsis_irradiance[40:length(t_csim)]))[1]
    csim_scale[i] <- (tsis_irradiance / csim_irradiance)[40:length(t_csim)][j] # around 6/2019
  }
}

### -------- binned ssi comparison ----------
{
  for (i in 1:4) {
    df <- data.frame(
      t = t_csim,
      tsis = binned_tsis[i,t_csim],
      siap = binned_siap[i,t_csim],
      siap_low = binned_siap[i,t_csim] - binned_siap_q[i,t_csim],
      siap_high = binned_siap[i,t_csim] + binned_siap_q[i,t_csim],
      gp = binned_gp[i,t_csim],
      gp_low = binned_gp_intervals[1,i,t_csim],
      gp_high = binned_gp_intervals[2,i,t_csim],
      marss = binned_marss[i,t_csim],
      marss_low = binned_marss[i,t_csim] - binned_marss_q[i,t_csim],
      marss_high = binned_marss[i,t_csim] + binned_marss_q[i,t_csim],
      csim = binned_csim[i,t_csim] * csim_scale[i]
    )
    
    df <- df %>%
      pivot_longer(
        cols = c("csim", "siap", "gp", "marss","tsis"), 
        names_to = "method", 
        values_to = "y"
      ) %>%
      mutate(
        low = if_else(method == "gp", gp_low, 
                      if_else(method == "marss", marss_low, 
                              if_else(method == "siap", siap_low, NA))), 
        high = if_else(method == "gp", gp_high, 
                       if_else(method == "marss", marss_high, 
                               if_else(method == "siap", siap_high, NA)))  
      )
    labels <- c(csim = paste0("CSIM*", round(csim_scale[i], 3)), siap = "SIAP",
                gp = "GP", marss = "MARSS", tsis = "TSIS")
    colors <- c(csim = "black", tsis = "grey", gp = "#619CFF", marss = "#B79F00", siap = "#F8766D")
    p <- ggplot(df, aes(x = t, y = y, group = method)) +
      geom_line(aes(color = method)) +  # Line for mean values
      # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y=0
      geom_ribbon(data = subset(df, method != "tsis" & method != "csim"), aes(ymin = low, ymax = high, fill = method), alpha = 0.2) +
      labs(title = c("210-300nm", "300-400nm", "700-1000nm", "1000-1300nm")[i],
           x = "Date (YYYY-MM-DD)", y = expression(Irradiance ~ "(" * mW ~ m^{-2} * ")"), color = "Method", fill = "Uncertainties") +
      scale_color_manual(values = colors, labels = labels) +
      scale_fill_manual(values = colors, labels = labels) +
      scale_x_discrete(breaks = c("2019-03-27", "2019-05-27", "2019-07-27", "2019-09-27", "2019-11-27")) +
      theme_minimal() +
      guides(fill = "none")
    ggsave(paste0("./output/realdata/binned_traintest_",i,"_.jpeg"), plot = p, width = 6, height = 3, dpi = 300)
  }
}

### -------- binned ssi comparison (zoomed in) ----------
{
  for (i in 1:4) {
    highlight_points <- t[which(apply(siap.res$x, 2, function(v) all(is.na(v))))]
    df <- data.frame(
      t = t_csim,
      tsis = binned_tsis[i,t_csim],
      siap = binned_siap[i,t_csim],
      siap_low = binned_siap[i,t_csim] - binned_siap_q[i,t_csim],
      siap_high = binned_siap[i,t_csim] + binned_siap_q[i,t_csim],
      gp = binned_gp[i,t_csim],
      gp_low = binned_gp_intervals[1,i,t_csim],
      gp_high = binned_gp_intervals[2,i,t_csim],
      marss = binned_marss[i,t_csim],
      marss_low = binned_marss[i,t_csim] - binned_marss_q[i,t_csim],
      marss_high = binned_marss[i,t_csim] + binned_marss_q[i,t_csim],
      csim = binned_csim[i,t_csim] * csim_scale[i]
    )
    df <- df %>%
      pivot_longer(
        cols = c("csim", "siap", "gp", "marss","tsis"), 
        names_to = "method", 
        values_to = "y"
      ) %>%
      mutate(
        low = if_else(method == "gp", gp_low, 
                      if_else(method == "marss", marss_low, 
                              if_else(method == "siap", siap_low, NA))), 
        high = if_else(method == "gp", gp_high, 
                       if_else(method == "marss", marss_high, 
                               if_else(method == "siap", siap_high, NA)))  
      )
    labels <- c(csim = paste0("CSIM*", round(csim_scale[i], 3)), siap = "SIAP",
                gp = "GP", marss = "MARSS", tsis = "TSIS")
    colors <- c(csim = "black", tsis = "grey", gp = "#619CFF", marss = "#B79F00", siap = "#F8766D")
    p <- ggplot(df, aes(x = t, y = y, group = method)) +
      geom_line(aes(color = method)) +  # Line for mean values
      # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y=0
      geom_ribbon(data = subset(df, method != "tsis" & method != "csim"), aes(ymin = low, ymax = high, fill = method), alpha = 0.2) +
      labs(title = c("210-300nm", "300-400nm", "700-1000nm", "1000-1300nm")[i],
           x = "Date (YYYY-MM-DD)", y = expression(Irradiance ~ "(" * mW ~ m^{-2} * ")"), color = "Method", fill = "Uncertainties") +
      scale_color_manual(values = colors, labels = labels) +
      scale_fill_manual(values = colors, labels = labels) +
      scale_x_discrete(breaks = c("2019-03-27", "2019-05-27", "2019-07-27", "2019-09-27", "2019-11-27")) +
      theme_minimal() +
      guides(fill = "none") +
      geom_point(data = df %>% filter(t %in% highlight_points), 
                 aes(x = t, y = y, color = method), 
                 size = 1, stroke = 0.5)  # Increase size for emphasis
    ggsave(paste0("./output/realdata/binned_traintest_",i,"_large.jpeg"), plot = p, width = 12, height = 3, dpi = 300)
  }
}

## ---------- hyperparametr tuning result ----------
### alpha path
{
  err <- readRDS("./output/realdata/alpha_path.rds")
  alpha_ <- readRDS("./output/realdata/alpha.rds")
  fig <- ggplot(data = err, aes(x = alpha, y = err)) + 
    geom_line() +
    labs(x = expression(alpha), y = "CV downtime MRAE") +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal(base_size = 16) +
    geom_point(data = err %>% filter(alpha == alpha_), 
               aes(x = alpha, y = err), 
               size = 2)
  
  ggsave(paste0("./output/realdata/hyperparam_path_alpha.jpeg"), plot = fig, width = 10, height = 3, dpi = 300)
}

# ------------- output reconstruction data ----------
{
  write.csv(siap_imp1, file = "output/realdata/siap_ssi_reconstruced.csv")
  write.csv(siap_q, file = "output/realdata/siap_ssi_95_uq.csv")
}
