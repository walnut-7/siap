source("code/read_data_functions.R")
source("code/plot_functions.R")
source("code/preprocess_functions.R")
source("code/postprocess_functions.R")
source("code/figure_config.R")

library(tidyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(patchwork)

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
  
  write.table(rae_stats, file = "./output/realdata/rae_stats_traintest.csv", sep = "\t")

  long_stats <- rae_stats %>%
    pivot_longer(
      cols = starts_with("siap") | starts_with("gp") | starts_with("marss") | starts_with("si"),
      names_to = c("method", ".value"),
      names_pattern = "(.*)_(.*)"
    )
  type_labels <- c(w = "downtime", o = "scattered")

  p <- ggplot(long_stats, aes(x = method, y = mean, col = method, shape = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +  # Scatter points
    geom_errorbar(aes(ymin = mean - q, ymax = mean + q), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      x = "Missingness ratio",
      y = "relative MRAE margin w.r.t. row-wise mean",
      color = "Method"
    ) +
    facet_wrap(~type, nrow = 1, labeller = as_labeller(type_labels)) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_shape_manual(values = method_shapes, labels = method_labels, name = "Method") +
    theme_minimal(base_size = 12) +
    guides(color = guide_legend(byrow = TRUE, nrow = 1, title.position = "left")) +
    theme(legend.position = "top") +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      axis.title.y.right = element_blank()
    )

  ggsave(paste0("./output/realdata/mrae_margin_to_mean_imputation.jpeg"), plot = p, width = 6, height = 4, dpi = 300)
}

## ----------- runtime -------------
# {
#   siap.res$siap$fit1$time + siap.res$siap$fit2$time
#   gp.res$runtime$runtime.gp + gp.res$runtime$runtime.pred + gp.res$runtime$runtime.sims
#   as.numeric(sum(marss.res$runtime$runtime.marss) + sum(marss.res$runtime$runtime.boot), units = "mins")
#   as.numeric(si_runtime, units = "mins")
# }

## ---------- UQ --------------
{
  coverage_df <- data.frame(siap = c(siap.res$ave_coverage, siap.res$ave_coverage.w, siap.res$ave_coverage.o),
                            gp = c(gp.res$ave_coverage, gp.res$ave_coverage.w, gp.res$ave_coverage.o),
                            marss = c((marss.res$ttl_coverage.w + marss.res$ttl_coverage.o) / (marss.res$n.w + marss.res$n.o),
                                      marss.res$ttl_coverage.w/marss.res$n.w,
                                      marss.res$ttl_coverage.o/marss.res$n.o))
  
  rownames(coverage_df) <- c("overall", "downtime", "scattered")
  write.table(coverage_df, file = "./output/realdata/coverage_traintest.csv", sep = "\t")
}

## --------- band-integrated ssi ------------
### helpers function
{
  integrate_ssi <- function(ssi) {
    # - ssi: Matrix. Row names must be the wavelength
    wvl <- as.numeric(rownames(ssi))
    area <- colSums((ssi[-1,,drop=F] + ssi[-nrow(ssi),,drop=F]) * (wvl[-1] - wvl[-length(wvl)]) / 2)
    area <- matrix(area, 1)
    rownames(area) <- paste0(round(range(wvl),1), collapse = "~")
    colnames(area) <- colnames(ssi)
    return(area)
  }
  
  binned <- function(ssi) {
    wvl <- as.numeric(rownames(ssi))
    rbind(
      integrate_ssi(ssi[which(wvl>210)[1]:which(wvl>300)[1], , drop=F]) / 1000,
      integrate_ssi(ssi[which(wvl>300)[1]:which(wvl>400)[1], , drop=F]) / 1000,
      integrate_ssi(ssi[which(wvl>400)[1]:which(wvl>700)[1], , drop=F]) / 1000,
      integrate_ssi(ssi[which(wvl>700)[1]:which(wvl>1000)[1], , drop=F]) / 1000,
      integrate_ssi(ssi[which(wvl>1000)[1]:which(wvl>=2399)[1], , drop=F]) / 1000
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

  binned_gp_sims <- array(dim=c(5, ncol(x), dim(gp_sims)[3]))
  for (i in 1:dim(gp_sims)[3]) {
    binned_gp_sims[,,i] <- binned(gp_sims[,,i])
  }
  binned_gp_intervals <- apply(binned_gp_sims, c(1,2), quantile, probs = c(0.025,0.975), na.rm = T)
  dimnames(binned_gp_intervals) <- list(
    rows = c("low", "high"),
    cols = c("210~300.1", "300.1~400.2", "400.2~701.1", "701.1~1002.9", "1002.9~2399.0"),
    slices = colnames(x)
  )

  gp_q <- gp.res$x_sd * qnorm(0.975)
  marss_q <- marss.res$x_sd * qnorm(0.975)

  rownames(gp_q) <- wvl
  binned_gp_q <- binned(gp_q)

  binned_marss <- binned(marss_imp)

  rownames(marss_q) <- wvl
  binned_marss_q <- binned(marss_q) 

  # --------- tsis
  binned_tsis <- binned(x)

  # -------- csim
  binned_csim <- matrix(NA,5,ncol(x))
  rownames(binned_csim) <- c("210~300.1", "300.1~400.2", "400.2~701.1", "701.1~1002.9", "1002.9~2399.0")
  colnames(binned_csim) <- colnames(x)
  binned_csim[,t_csim] <- binned(csim_uniq)

  # ---------- baseline
  pre_par <- preprocess_func(x, box.cox = F)
  x_ready <- pre_par$x_ready
  x_bsl <- x_ready
  x_bsl[which(is.na(x_ready), arr.ind = T)] <- pre_par$nmlz_mu[which(is.na(x_ready), arr.ind = T)[, 1]]
  x_bsl <- postprocess_func(x_bsl, pre_par = pre_par, x)
  rownames(x_bsl) <- wvl
  colnames(x_bsl) <- colnames(x)
  binned_baseline <- binned(x_bsl)

}

### scale csim
{
  wvl_ranges_ls <- list(c(210,300), c(300,400), c(400, 700), c(700,1000), c(1000,2400))
  csim_scale <- rep(NA, 5)

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
  highlight_points <- t[which(apply(siap.res$x, 2, function(v) all(is.na(v))))]
  full_t <- # t %>% as.character()
    # seq.Date(min(ymd(t_csim)), max(ymd(t_csim)), by = "day") %>% as.character()
    seq.Date(t[1], max(ymd(t_csim)), by = "day") %>% as.character()

  align_to_full_t <- function(values) {
    out <- rep(NA_real_, length(full_t))
    out[match(t_csim, full_t)] <- as.numeric(values)
    out
  }

  make_binned_panel <- function(i, show_missing_points = FALSE) {
    df <- data.frame(
      t = full_t %>% as.Date(),
      tsis = binned_tsis[i, full_t],
      tsis_baseline = binned_baseline[i, full_t],
      siap = binned_siap[i, full_t],
      siap_low = binned_siap[i, full_t] - binned_siap_q[i, full_t],
      siap_high = binned_siap[i, full_t] + binned_siap_q[i, full_t],
      gp = binned_gp[i, full_t],
      gp_low = binned_gp_intervals[1, i, full_t],
      gp_high = binned_gp_intervals[2, i, full_t],
      marss = binned_marss[i, full_t],
      marss_low = binned_marss[i, full_t] - binned_marss_q[i, full_t],
      marss_high = binned_marss[i, full_t] + binned_marss_q[i, full_t],
      csim = binned_csim[i, full_t] * csim_scale[i]
    )
    
    df <- df %>%
      pivot_longer(
        cols = c("csim", "siap", "gp", "marss", "tsis", "tsis_baseline"),
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
    
    realdata_labels <- c(csim = paste0("CSIM*", round(csim_scale[i], 3)), siap = "SIAP",
            gp = "GP", marss = "MARSS", tsis = "TSIS", tsis_baseline = "TSIS-baseline")

    p <- ggplot(df, aes(x = t, y = y, group = method)) +
      geom_line(aes(color = method)) +
      geom_ribbon(
        data = subset(df, method %in% c("gp", "marss", "siap")),
        aes(ymin = low, ymax = high, fill = method),
        alpha = 0.2
      ) +
      labs(
        title = c("210-300nm", "300-400nm", "400-700nm", "700-1000nm", "1000-2400nm")[i],
        x = "Date (YYYY-MM-DD)",
        y = expression(Irradiance ~ "(" * W ~ m^{-2} * ")"),
        color = "Method",
        fill = "Uncertainties"
      ) +
      scale_color_manual(values = realdata_colors, labels = realdata_labels) +
      scale_fill_manual(values = realdata_colors, labels = realdata_labels) +
      scale_x_date(
        breaks = as.Date(c("2019-03-27", "2019-05-27", "2019-07-27", "2019-09-27", "2019-11-27")),
        date_labels = "%Y-%m-%d"
      ) +
      theme_minimal() +
      guides(fill = "none")
    
    if (show_missing_points) {
      p <- p + geom_point(
        data = df %>% filter(t %in% highlight_points),
        aes(x = t, y = y, color = method),
        size = 1,
        stroke = 0.5
      )
    }
    
    p
  }


  for (i in 1:5) {
    p <- make_binned_panel(i, show_missing_points = FALSE)
    ggsave(paste0("./output/realdata/binned_traintest_",i,"_.jpeg"), plot = p, width = 6, height = 3, dpi = 300)
  }
}

### -------- binned ssi comparison (combined) ----------
{
  highlight_points <- t[which(apply(siap.res$x, 2, function(v) all(is.na(v))))]
  full_t <- seq.Date(t[1], t[500], by = "day") %>% as.character()

  align_to_full_t <- function(values) {
    out <- rep(NA_real_, length(full_t))
    out[match(t_csim, full_t)] <- as.numeric(values)
    out
  }

  make_binned_panel <- function(i, show_missing_points = FALSE) {
    df <- data.frame(
      t = full_t %>% as.Date(),
      tsis = binned_tsis[i, full_t],
      siap = binned_siap[i, full_t],
      siap_low = binned_siap[i, full_t] - binned_siap_q[i, full_t],
      siap_high = binned_siap[i, full_t] + binned_siap_q[i, full_t],
      gp = binned_gp[i, full_t],
      gp_low = binned_gp_intervals[1, i, full_t],
      gp_high = binned_gp_intervals[2, i, full_t],
      marss = binned_marss[i, full_t],
      marss_low = binned_marss[i, full_t] - binned_marss_q[i, full_t],
      marss_high = binned_marss[i, full_t] + binned_marss_q[i, full_t],
      csim = binned_csim[i, full_t] * csim_scale[i]
    )
    
    df <- df %>%
      pivot_longer(
        cols = c("csim", "siap", "gp", "marss", "tsis"),
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
    realdata_labels <- c(csim = paste0("CSIM*", round(csim_scale[i], 3)), siap = "SIAP",
            gp = "GP", marss = "MARSS", tsis = "TSIS", tsis_baseline = "TSIS-baseline")
    p <- ggplot(df, aes(x = t, y = y, group = method)) +
      geom_ribbon(
        data = subset(df, method != "tsis" & method != "csim"),
        aes(ymin = low, ymax = high, fill = method),
        alpha = 0.2
      ) +
      geom_line(data = subset(df, method != "csim"), aes(color = method)) +
      geom_line(data = subset(df, method == "csim"), aes(color = method)) +
      labs(
        title = c("210-300nm", "300-400nm", "400-700nm", "700-1000nm", "1000-2400nm")[i],
        x = "Date (YYYY-MM-DD)",
        y = expression(Irradiance ~ "(" * W ~ m^{-2} * ")"),
        color = "Method",
        fill = "Uncertainties"
      ) +
      scale_color_manual(values = realdata_colors, labels = realdata_labels) +
      scale_fill_manual(values = realdata_colors, labels = realdata_labels) +
      scale_x_date(
        breaks = as.Date(c("2018-03-27", "2018-09-27", "2019-03-27", "2019-09-27")),
        date_labels = "%Y-%m-%d"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 9, face = "bold")) +
      guides(fill = "none")

    if (show_missing_points) {
      pts <- df %>% filter(t %in% highlight_points)
      p <- p +
        geom_point(
          data = subset(pts, method != "csim"),
          aes(x = t, y = y, color = method),
          size = 1,
          stroke = 0.5
        ) +
        geom_point(
          data = subset(pts, method == "csim"),
          aes(x = t, y = y, color = method),
          size = 1,
          stroke = 0.5
        )
    }

    p
  }

  p1 <- make_binned_panel(1, show_missing_points = FALSE)
  p3 <- make_binned_panel(3, show_missing_points = FALSE)
  p4 <- make_binned_panel(4, show_missing_points = FALSE)
  p5 <- make_binned_panel(5, show_missing_points = FALSE)
  p2 <- make_binned_panel(2, show_missing_points = TRUE)

  p1 <- p1 + labs(x = NULL)
  p3 <- p3 + labs(x = NULL, y = NULL)
  p5 <- p5 + labs(x = NULL, y = NULL)
  p4 <- p4 + labs(x = NULL)

  p_combined <- ((p1 + theme(legend.position = "none")) + (p3 + theme(legend.position = "none"))) /
    ((p4 + theme(legend.position = "none")) + (p5 + theme(legend.position = "none"))) /
    (p2 + theme(legend.position = "bottom"))

  ggsave("./output/realdata/binned_traintest_combined.jpeg", plot = p_combined, width = 8, height = 6, dpi = 300)
}

### ------- gp raw ------
{
  x_tr <- gp.res$x
  pre_par_tr <- preprocess_func(x_tr, box.cox = F)
  x_df <- data.frame(
    row = rep(1:nrow(pre_par_tr$x_ready), times = ncol(pre_par_tr$x_ready)),
    col = rep(1:ncol(pre_par_tr$x_ready), each = nrow(pre_par_tr$x_ready)),
    value = as.vector(pre_par_tr$x_ready)
  )
  is_missing <- is.na(x_df$value)
  obs_and_pred <- x_df$value
  obs_and_pred[is_missing] <- gp.res$pred
  x_pred_raw <- array(obs_and_pred, c(nrow(pre_par_tr$x_ready), ncol(pre_par_tr$x_ready)))
  rownames(x_pred_raw) <- wvl
  binned_raw <- binned(x_pred_raw)
  colnames(binned_raw) <- colnames(x_ready)

  locs <- cbind(x_df$row, x_df$col)
  z <- matrix(rep(1, d2), ncol = 1) # n x 1
  exog <- z %*% gp.res$gp$betahat
  exog <- array(exog, c(nrow(pre_par_tr$x_ready), ncol(pre_par_tr$x_ready)))
  rownames(exog) <- wvl
  colnames(exog) <- colnames(x_ready)
  binned_exog <- binned(exog)

  pre_par <- preprocess_func(x, box.cox = F)
  x_ready <- pre_par$x_ready
  binned_tsis_preprocess <- binned(x_ready)

  i <- 2

  t_grid <- # colnames(x_ready)
    seq.Date(min(ymd(t_csim)), max(ymd(t_csim)), by = "day") %>% as.character()
  df <- data.frame(
    t = t_grid %>% as.Date(),
    tsis_preprocess = binned_tsis_preprocess[i, t_grid],
    exog = binned_exog[i, t_grid],
    gp_raw = binned_raw[i, t_grid]
    )

  df <- df %>%
    pivot_longer(
      cols = c("tsis_preprocess", "gp_raw", "exog"),
      names_to = "method",
      values_to = "y"
    ) %>%
    mutate(method = recode(method,
      tsis_preprocess = "tsis",
      gp_raw = "gp"
    ))

  plot_colors <- c(realdata_colors[c("tsis", "gp")], exog = "darkorange")
  plot_labels <- c(gp = "GP", tsis = "TSIS", exog = "Mean estimate of GP")
  band_titles <- c("210-300nm", "300-400nm", "400-700nm", "700-1000nm", "1000-2400nm")

  highlight_points <- colnames(pre_par_tr$x_ready)[apply(pre_par_tr$x_ready, 2, function(v) all(is.na(v)))]
  highlight_points <- intersect(highlight_points, t_grid)
  df_missing <- data.frame(
    t = as.Date(highlight_points),
    y = binned_raw[i, highlight_points],
    method = "gp"
  )

  p <- ggplot(df, aes(x = t, y = y, color = method)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(data = df_missing, aes(x = t, y = y, color = method),
              shape = 16, size = 1.5, inherit.aes = FALSE) +
    scale_color_manual(values = plot_colors, labels = plot_labels) +
    scale_x_date(
      breaks = as.Date(c("2019-03-27", "2019-05-27", "2019-07-27", "2019-09-27", "2019-11-27")),
      date_labels = "%b %Y"
    ) +
    labs(
      title = band_titles[i],
      x = "Date",
      y = "Preprocessed irradiance",
        # expression(Delta ~ Irradiance ~ "(standardized)"),
      color = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(paste0("./output/realdata/binned_gp_raw_",i , ".jpeg"), plot = p, width = 8, height = 4, dpi = 300)
}

### -------- downtime missingess type summary -------
{
  t_grid <- colnames(x_ready)
  binned_tsis_preprocess <- binned(x_ready)

  hp_dates <- sort(as.Date(highlight_points))
  t_dates <- as.Date(t_grid)

  gaps <- c(TRUE, diff(hp_dates) > 1)
  chunks_tbl <- data.frame(date = hp_dates, chunk = cumsum(gaps)) %>%
    group_by(chunk) %>%
    summarise(start = min(date), end = max(date), .groups = "drop") %>%
    mutate(width = as.integer(end - start) + 1L)

  # for each band, classify each chunk by the sign of its flanking observed values
  chunk_summary <- lapply(1:5, function(i) {
    classify_chunk <- function(start, end) {
      left_t  <- as.character(t_dates[t_dates < start])
      right_t <- as.character(t_dates[t_dates > end])
      left_cands  <- left_t[!is.na(binned_tsis_preprocess[i, left_t])]
      right_cands <- right_t[!is.na(binned_tsis_preprocess[i, right_t])]
      if (length(left_cands) == 0 || length(right_cands) == 0) return(NA_character_)
      lv <- binned_tsis_preprocess[i, left_cands[length(left_cands)]]
      rv <- binned_tsis_preprocess[i, right_cands[1]]
      if (is.na(lv) || is.na(rv)) return(NA_character_)
      if (lv > 0 && rv > 0) "pos-pos"
      else if (lv < 0 && rv < 0) "neg-neg"
      else "neg-pos"
    }
    types <- mapply(classify_chunk, chunks_tbl$start, chunks_tbl$end)
    widths <- chunks_tbl$width
    n <- length(types)
    valid <- !is.na(types)
    w_tot <- sum(widths[valid])
    data.frame(
      i = i,
      n_chunks = n,
      prop_pos_pos = sum(types == "pos-pos", na.rm = TRUE) / n,
      prop_neg_neg = sum(types == "neg-neg", na.rm = TRUE) / n,
      prop_neg_pos = sum(types == "neg-pos", na.rm = TRUE) / n,
      wprop_pos_pos = sum(widths[types == "pos-pos" & valid]) / w_tot,
      wprop_neg_neg = sum(widths[types == "neg-neg" & valid]) / w_tot,
      wprop_neg_pos = sum(widths[types == "neg-pos" & valid]) / w_tot
    )
  }) %>% bind_rows()

  write.csv(chunk_summary, file = "output/realdata/chunk_summary.csv")
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
