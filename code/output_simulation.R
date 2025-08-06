source("code/read_data_functions.R")

dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt

library(dplyr)
library(ggplot2)
library(ggbreak)
library(scales)
library(ggpubr)

# ------------ Figures -------------
method_labels <- c(
  siap = "SIAP",
  si = "SI",
  si_trend = "SI detrended",
  si_trend_cov = "SI detrended w/ cov (S1)",
  sia = "SS SIA",
  sia_trend = "SS SIA detrended (S2)",
  s1_si_trend = "S1 + SS SI detrended",
  gp = "GP",
  marss = "MARSS"
)

## -------- rel mrae margin, grouped by p_mis --------
for (pdt in c(0.1, 0.3, 0.5)) {
  df <- readRDS(file = paste0("./output/simulation/overall_rel_mrae_margin_", 10*pdt, ".rds"))
  df_stats <- df %>%
    group_by(type, method) %>%
    summarise(
      mean_mrae = mean(rel_margin_mrae, na.rm = T),
      wid_mrae = sd(rel_margin_mrae, na.rm = T) / sqrt(sum(!is.na(rel_margin_mrae))) * qnorm(0.975),
      .groups = "drop"
    )
  p <- ggplot(df_stats[!df_stats$method%in%c("sia","si", "sia_trend","gp", "marss"),], aes(x = type, y = mean_mrae, col = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 1) +  # Scatter points
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae),
                  width = 0.1, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = paste0(pdt*10, "0% downtime and scattered missingness ratio"),
      x = "Type",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    theme_minimal() +
    scale_color_discrete(labels = method_labels) +
    scale_x_discrete(labels = c("o" = "scattered missingness", "w" = "downtime missingness")) #+
  #ylim(c(-0.01,0.9))
  ggsave(paste0("./output/simulation/rel_mrae_margin_to_spline_positive_",10*pdt,".jpeg"), plot = p, width = 6, height = 6, dpi = 300)
  
  p <- ggplot(df_stats[df_stats$method%in%c("sia","si", "sia_trend","gp", "marss"),], aes(x = type, y = mean_mrae, col = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 1) +  # Scatter points
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae),
                  width = 0.1, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = paste0(pdt*10, "0% downtime and scattered missingness ratio"),
      x = "Type",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    theme_minimal() +
    scale_y_break(c(-110, -1050)) +
    scale_color_discrete(labels = method_labels) +
    scale_x_discrete(labels = c("o" = "scattered missingness", "w" = "downtime missingness")) #+
  #ylim(c(-1160,0.42))
  ggsave(paste0("./output/simulation/rel_mrae_margin_to_spline_negative_",10*pdt,".jpeg"), plot = p, width = 6, height = 6, dpi = 300)
}

## -------- rel mrae margin, grouped by missingness type ---------
{
  df1 <- readRDS(file = paste0("./output/simulation/overall_rel_mrae_margin_1.rds"))
  df2 <- readRDS(file = paste0("./output/simulation/overall_rel_mrae_margin_3.rds"))
  df3 <- readRDS(file = paste0("./output/simulation/overall_rel_mrae_margin_5.rds"))
  df1$pdt = 0.1
  df2$pdt = 0.3
  df3$pdt = 0.5
  df = rbind(df1, df2, df3)
  df$pdt <- as.factor(df$pdt)
  
  df_stats <- df %>%
    group_by(type, method, pdt) %>%
    summarise(
      mean_mrae = mean(rel_margin_mrae, na.rm = T),
      wid_mrae = sd(rel_margin_mrae, na.rm = T) / sqrt(sum(!is.na(rel_margin_mrae))) * qnorm(0.975),
      .groups = "drop"
    )
  
  p <- ggplot(df_stats[!df_stats$method%in%c("sia","si", "sia_trend","gp","marss")&df_stats$type=="w",], aes(x = pdt, y = mean_mrae, col = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.5) +  # Scatter points
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = "downtime",
      x = "Missingness ratio",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    scale_color_discrete(labels = method_labels) +
    theme_minimal(base_size = 14) #+
  #ylim(c(-0.01,0.03))
  
  ggsave(paste0("./output/simulation/rel_mrae_margin_to_spline_positive_w.jpeg"), plot = p, width = 6, height = 3, dpi = 300)
  
  p <- ggplot(df_stats[!df_stats$method%in%c("sia","si", "sia_trend","gp","marss")&df_stats$type=="o",], aes(x = pdt, y = mean_mrae, col = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.5) +  # Scatter points
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = "scattered",
      x = "Missingness ratio",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    scale_color_discrete(labels = method_labels) +
    #ylim(c(0.84,0.9)) +
    theme_minimal(base_size = 14)
  
  ggsave(paste0("./output/simulation/rel_mrae_margin_to_spline_positive_o.jpeg"), plot = p, width = 6, height = 3, dpi = 300)
  
  p <- ggplot(df_stats[df_stats$method%in%c("sia","si", "sia_trend","gp","marss")&df_stats$type=="w",], aes(x = pdt, y = mean_mrae, col = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.5) +  # Scatter points
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = "downtime",
      x = "Missingness ratio",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    theme_minimal(base_size = 14) +
    scale_y_break(c(-50, -1050)) +
    scale_color_discrete(labels = method_labels) +
    ylim(c(-1100,0.42))
  
  ggsave(paste0("./output/simulation/rel_mrae_margin_to_spline_negative_w.jpeg"), plot = p, width = 6, height = 3, dpi = 300)
  
  p <- ggplot(df_stats[df_stats$method%in%c("sia","si", "sia_trend","gp","marss")&df_stats$type=="o",], aes(x = pdt, y = mean_mrae, col = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 0.5) +  # Scatter points
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = "scattered",
      x = "Missingness ratio",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    scale_color_discrete(labels = method_labels) +
    theme_minimal(base_size = 14) #+
  #ylim(c(-150,0.42))
  
  ggsave(paste0("./output/simulation/rel_mrae_margin_to_spline_negative_o.jpeg"), plot = p, width = 6, height = 3, dpi = 300)
}


## ------------- mrae v.s. wvl -----------------
for (pdt in c(0.1, 0.3, 0.5)) {
  mrae_stats <- readRDS(file = paste0("./output/simulation/mrae_stats_",pdt*10,".rds"))
  mrae_stats$nrepl <- 100
  if (pdt == 0.5) {
    mrae_stats$nrepl[which(mrae_stats$method=="marss")] <- 85
  }
  plot_data <- mrae_stats %>%
    mutate(
      lower_w = mean_mrae.w - 1.96 * sd_mrae.w / sqrt(nrepl),
      upper_w = mean_mrae.w + 1.96 * sd_mrae.w / sqrt(nrepl),
      lower_o = mean_mrae.o - 1.96 * sd_mrae.o / sqrt(nrepl),
      upper_o = mean_mrae.o + 1.96 * sd_mrae.o / sqrt(nrepl)
    )
  p1 <- ggplot(plot_data) +
    geom_ribbon(aes(x = wvl, ymin = lower_w, ymax = upper_w, fill = method), alpha = 0.2) +
    geom_line(aes(x = wvl, y = mean_mrae.w, color = method, group = method)) +
    labs(
      x = "Wavelength (nm)",
      y = "Mean Relative Absolute Error",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_discrete(labels = method_labels) +
    scale_fill_discrete(labels = method_labels) +
    scale_y_break(c(0.1, 0.8), scales = 0.1) # Define the gap in the y-axis
  
  ggsave(paste0("./output/simulation/mrae_w_",10*pdt,"_.jpeg"), plot = p1, width = 6, height = 3, dpi = 300)
  
  p2 <- ggplot(plot_data) +
    geom_ribbon(aes(x = wvl, ymin = lower_o, ymax = upper_o, fill = method), alpha = 0.2) +
    geom_line(aes(x = wvl, y = mean_mrae.o, color = method, group = method)) +
    labs(
      x = "Wavelength (nm)",
      y = "Mean Relative Absolute Error",
      color = "Method",
      fill = "Method"
    ) +
    scale_color_discrete(labels = method_labels) +
    scale_fill_discrete(labels = method_labels) +
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10() 
  
  ggsave(paste0("./output/simulation/mrae_o_",10*pdt,"_.jpeg"), plot = p2, width = 6, height = 3, dpi = 300)
}

## ------------- interval length v.s. wvl ----------
for (pdt in c(0.1, 0.3, 0.5)) {
  suffix <- c("", "ver1.0_")[1]
  wvl <- as.numeric(rownames(xt))
  length_stats <- readRDS(file = paste0("./output/simulation/length_stats_", suffix, pdt*10,".rds"))
  ex <- apply(xt, 1, mean)
  
  plot_data <- length_stats 
  
  p1 <- ggplot(plot_data) +
    #geom_ribbon(aes(x = wvl, ymin = lower_w, ymax = upper_w, fill = method), alpha = 0.2) +
    geom_line(aes(x = wvl, y = mean_length.w/ex, color = method, group = method)) + # half interval length (i.e. one-sided)
    labs(
      x = "Wavelength (nm)",
      y = "Relative Interval Width",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10(breaks = scales::log_breaks(n = 10)) +
    ggbreak::scale_y_break(c(0.1, .8)) # Define the gap in the y-axis
  ggsave(paste0("./output/simulation/rel_length_w_", suffix, 10*pdt,"_.jpeg"), plot = p1, width = 6, height = 3, dpi = 300)
  
  p1 <- ggplot(plot_data) +
    geom_line(aes(x = wvl, y = mean_length.w, color = method, group = method)) +
    labs(
      x = "Wavelength (nm)",
      y = "Interval Width",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10(breaks = scales::log_breaks(n = 10))
  
  ggsave(paste0("./output/simulation/length_w_", suffix, 10*pdt,"_.jpeg"), plot = p1, width = 6, height = 3, dpi = 300)
  
  p2 <- ggplot(plot_data) +
    geom_line(aes(x = wvl, y = mean_length.o/ex, color = method, group = method)) +
    labs(
      x = "Wavelength (nm)",
      y = "Relative Interval Width",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10(breaks = scales::log_breaks(n = 10))
  
  ggsave(paste0("./output/simulation/rel_length_o_", suffix, 10*pdt,"_.jpeg"), plot = p2, width = 6, height = 3, dpi = 300)
  
  p2 <- ggplot(plot_data) +
    geom_line(aes(x = wvl, y = mean_length.o, color = method, group = method)) +
    labs(
      x = "Wavelength (nm)",
      y = "Interval Width",
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    scale_y_log10() +
    scale_x_log10(breaks = scales::log_breaks(n = 10))
  
  ggsave(paste0("./output/simulation/length_o_", suffix, 10*pdt,"_.jpeg"), plot = p2, width = 6, height = 3, dpi = 300)
}


## ------------- average coverage v.s. wvl --------------
{
  suffix <- c("", "ver1.0_")[1]
  df_combined <- NULL
  for (i in c(1,2,3)) {
    pdt <- c(0.1, 0.3, 0.5)[i]
    df <- readRDS(file = paste0("./output/simulation/coverage_wvl_stats_", suffix, 10*pdt, ".rds"))
    df$wvl <- as.numeric(df$wvl)
    
    df$nrepl <- 100
    if (pdt == 0.5) {
      df$nrepl[which(df$method=="marss")] <- 85
    }
    df <- df %>%
      mutate(
        lower_overall = overall - 1.96 * sd_overall / sqrt(nrepl),
        upper_overall = overall + 1.96 * sd_overall / sqrt(nrepl),
        lower_w = w - 1.96 * sd_w / sqrt(nrepl),
        upper_w = w + 1.96 * sd_w / sqrt(nrepl),
        lower_o = o - 1.96 * sd_o / sqrt(nrepl),
        upper_o = o + 1.96 * sd_o / sqrt(nrepl)
      )
    df_combined <- rbind(df_combined, mutate(df, pdt = pdt))
    p[[i]] <- ggplot(df, aes(x = wvl, y = overall, color = method, group = method)) +
      geom_ribbon(aes(x = wvl, ymin = lower_overall, ymax = upper_overall, fill = method), alpha = 0.2, color = NA) +
      geom_line() +  # Line for each method
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 1) + 
      labs(x = "Wavelength (nm)", y = "Mean coverage rate", title = "") +
      scale_x_log10(breaks = scales::log_breaks(n = 10)) + 
      scale_color_discrete(labels = method_labels) +
      scale_fill_discrete(labels = method_labels) +
      theme_minimal() +
      coord_cartesian(ylim = c(0, 1))
    
    ggsave(paste0("./output/simulation/coverage_wvl_", suffix, 10*pdt,".jpeg"), 
           plot = p[[i]] + scale_y_break(c(0.2, 0.6)), 
           width = 6, height = 3, dpi = 300)
    
    p1 <- ggplot(df, aes(x = wvl, y = o, color = method, group = method)) +
      geom_line() +  
      geom_ribbon(aes(x = wvl, ymin = lower_o, ymax = upper_o, fill = method), alpha = 0.2, color = NA) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 1) + 
      labs(x = "Wavelength (nm)", y = "Mean coverage rate", title = "") +
      scale_x_log10(breaks = scales::log_breaks(n = 10)) + 
      scale_color_discrete(labels = method_labels) +
      scale_fill_discrete(labels = method_labels) +
      theme_minimal() 
    
    ggsave(paste0("./output/simulation/coverage_wvl_o_", suffix, 10*pdt,".jpeg"), plot = p1, width = 6, height = 3, dpi = 300)
    
    p2 <- ggplot(df, aes(x = wvl, y = w, color = method, group = method)) +
      geom_line() +  
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", linewidth = 1) + 
      geom_ribbon(aes(x = wvl, ymin = lower_w, ymax = upper_w, fill = method), alpha = 0.2, color = NA) +
      labs(x = "Wavelength (nm)", y = "Mean coverage rate", title = "") +
      scale_x_log10(breaks = scales::log_breaks(n = 10)) + 
      scale_color_discrete(labels = method_labels) +
      scale_fill_discrete(labels = method_labels) +
      theme_minimal() 
    
    ggsave(paste0("./output/simulation/coverage_wvl_w_", suffix, 10*pdt,".jpeg"), plot = p1, width = 6, height = 3, dpi = 300)
  }
  combined_plot <- ggplot(df_combined, aes(x = wvl, y = overall, color = method, group = method)) +
    geom_ribbon(aes(ymin = lower_overall, ymax = upper_overall, fill = method), alpha = 0.2, color = NA) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    labs(x = "Wavelength (nm)", y = "Mean coverage rate", title = "") +
    scale_x_log10(breaks = scales::log_breaks(n = 8)) +
    scale_y_break(c(0.2, 0.6), scales = 2) +
    facet_wrap(~pdt, nrow = 1, 
               labeller = labeller(pdt = c("0.1" = "Missing Rate = 10%", 
                                           "0.3" = "Missing Rate = 30%", 
                                           "0.5" = "Missing Rate = 50%"))) +
    scale_color_discrete(labels = method_labels) +
    scale_fill_discrete(labels = method_labels) +
    theme_minimal(base_size = 14)
  
  p <- combined_plot +
    geom_curve(
      aes(x = 2100, y = 0.93, xend = 2300, yend = 0.9),
      arrow = arrow(length = unit(0.02, "npc")),
      curvature = -0.2,
      color = "#8A8E91"
    ) +
    annotate("text", x = 2300, y = 0.88, label = "SIAP-type Methods", hjust = 0.9, size = 4, color = "#8A8E91")
  
  ggsave("./output/simulation/coverage_wvl_annotated.jpeg", plot = p, width = 12, height = 5, dpi = 300)
  
  combined_plot <- ggplot(df_combined, aes(x = wvl, y = o, color = method, group = method)) +
    geom_ribbon(aes(ymin = lower_o, ymax = upper_o, fill = method), alpha = 0.2, color = NA) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    labs(x = "Wavelength (nm)", y = "Mean coverage rate", title = "") +
    scale_x_log10(breaks = scales::log_breaks(n = 8)) +
    scale_y_break(c(0.2, 0.6), scales = 2) +
    facet_wrap(~pdt, nrow = 1, 
               labeller = labeller(pdt = c("0.1" = "Missing Rate = 10%", 
                                           "0.3" = "Missing Rate = 30%", 
                                           "0.5" = "Missing Rate = 50%"))) +
    scale_color_discrete(labels = method_labels) +
    scale_fill_discrete(labels = method_labels) +
    theme_minimal(base_size = 14)
  
  ggsave("./output/simulation/coverage_wvl_o.jpeg", plot = combined_plot, width = 12, height = 5, dpi = 300)
  
  combined_plot <- ggplot(df_combined, aes(x = wvl, y = w, color = method, group = method)) +
    geom_ribbon(aes(ymin = lower_w, ymax = upper_w, fill = method), alpha = 0.2, color = NA) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
    labs(x = "Wavelength (nm)", y = "Mean coverage rate", title = "") +
    scale_x_log10(breaks = scales::log_breaks(n = 8)) +
    scale_y_break(c(0.2, 0.6), scales = 2) +
    facet_wrap(~pdt, nrow = 1, 
               labeller = labeller(pdt = c("0.1" = "Missing Rate = 10%", 
                                           "0.3" = "Missing Rate = 30%", 
                                           "0.5" = "Missing Rate = 50%"))) +
    scale_color_discrete(labels = method_labels) +
    scale_fill_discrete(labels = method_labels) +
    theme_minimal(base_size = 14)
  
  ggsave("./output/simulation/coverage_wvl_w.jpeg", plot = combined_plot, width = 12, height = 5, dpi = 300)
}

# ------------ Table -------------
## -------- average coverage rate -----------
{
  coverage_tb <-NULL
  
  for (pdt in c(0.1,0.3,0.5)) {
    coverage_stats <- readRDS(file = paste0("./output/simulation/coverage_stats_",10*pdt, ".rds"))
    coverage_stats$nrepl <- 100
    if (pdt == 0.5) {
      coverage_stats$nrepl[which(coverage_stats$method=="marss")] <- 85
    }
    coverage_stats = coverage_stats %>% 
      mutate(
        sd_coverage = sd_coverage / sqrt(nrepl),
        sd_coverage.w = sd_coverage.w / sqrt(nrepl),
        sd_coverage.o = sd_coverage.o / sqrt(nrepl)
      )
    coverage_stats$pdt = pdt
    coverage_tb <- rbind(coverage_tb, coverage_stats)
  }
  
  colnames(coverage_tb) <- c("mean_coverage", "se_coverage", "mean_coverage.w", "se_coverage.w", "mean_coverage.o", "se_coverage.o" , "method", "nrepl", "pdt")
  
  coverage_summary <- coverage_tb %>%
    group_by(pdt, method) %>%
    summarize(
      mean_coverage = round(mean_coverage, 4),
      se_coverage = round(se_coverage, 4),
      mean_coverage_w = round(mean_coverage.w, 4),
      se_coverage_w = round(se_coverage.w, 4),
      mean_coverage_o = round(mean_coverage.o, 4),
      se_coverage_o = round(se_coverage.o, 4)
    )
  
  write.table(coverage_summary, file = "./output/simulation/coverage_table.csv", row.names = F, sep = "\t")
}
