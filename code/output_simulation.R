source("code/read_data_functions.R")
source("./code/result_functions.R")
source('./code/figure_config.R')
library(dplyr)
library(ggplot2)
library(ggbreak)
library(tidyr)
library(ggh4x)

dl <- load_ssi(option = "synthetic", missing.option = "0")
xt <- dl$xt

# ----------- helper functions -------------
get_df_overall_mrae <- function(gp_diff_flag = T, marss_diff_flag = T){
  df1 <- readRDS(file = result_path_simu("overall_rel_mrae_margin", 0.1, gp_diff_flag, marss_diff_flag))
  df2 <- readRDS(file = result_path_simu("overall_rel_mrae_margin", 0.3, gp_diff_flag, marss_diff_flag))
  df3 <- readRDS(file = result_path_simu("overall_rel_mrae_margin", 0.5, gp_diff_flag, marss_diff_flag))
  
  df_trmf <- read.table(file = "./output/test9/overall_rel_mrae_margin_trmf.csv", sep = ',', header = T) %>%
    filter(type %in% c('w', 'o'))
  df_trmf$pdt <- as.factor(df_trmf$pdt)
  df_latc <- read.table(file = "./output/test9/overall_rel_mrae_margin_latc.csv", sep = ',', header = T) %>%
    filter(type %in% c('w', 'o'))
  df_latc$pdt <- as.factor(df_latc$pdt)
  
  df1$pdt = 0.1
  df2$pdt = 0.3
  df3$pdt = 0.5
  df = rbind(df1, df2, df3)
  df$pdt <- as.factor(df$pdt)
  
  df0 <- df %>% # baseline model mrae
    filter(method == 'siap') %>%
    mutate(mrae0 = mrae / (1 - rel_margin_mrae),
           repl = rep(rep(1:100, each = 2), 3)) %>%
    select(type, pdt, repl, mrae0)
  
  df_trmf <- df0 %>% 
    left_join(df_trmf, by = c('type', 'pdt', 'repl')) %>%
    mutate(rel_margin_mrae = (mrae0 - mrae) / mrae0) %>%
    arrange(pdt, repl, desc(type)) %>%
    select(type, mrae, method, rel_margin_mrae, pdt) 
  
  df_latc <- df0 %>% 
    inner_join(df_latc, by = c('type', 'pdt', 'repl')) %>%
    mutate(rel_margin_mrae = (mrae0 - mrae) / mrae0) %>%
    arrange(pdt, repl, desc(type)) %>%
    select(type, mrae, method, rel_margin_mrae, pdt) 
  
  df <- rbind(df, df_trmf, df_latc)

  return(df)
}


# ------------ Figures -------------
## ---------- mrae -------------
for (diff_flag in c(T,F)) {
  gp_diff_flag <- diff_flag
  marss_diff_flag <- diff_flag
  df <- get_df_overall_mrae(gp_diff_flag, marss_diff_flag)
  
  df_stats <- df %>%
    group_by(type, method, pdt) %>%
    summarise(
      mean_mrae = mean(rel_margin_mrae, na.rm = T),
      wid_mrae = sd(rel_margin_mrae, na.rm = T) / sqrt(sum(!is.na(rel_margin_mrae))) * qnorm(0.975),
      .groups = "drop"
    )
  
  p <- ggplot(df_stats[df_stats$method %in% subset_methods & df_stats$type=="w",],
                    aes(x = pdt, y = mean_mrae, col = method, shape = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +  # Scatter points
    geom_text(
      data = df_stats[df_stats$method %in% subset_methods & df_stats$type=="w" & df_stats$method=="siap",],
      aes(label = sprintf("%.3f", mean_mrae)),
      vjust = -.5, hjust = -.5,
      size = 3,
      show.legend = FALSE
    ) +
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = "downtime",
      x = "Missingness ratio",
      y = "relative MRAE margin",
      color = "Method"
    ) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_shape_manual(values = method_shapes, labels = method_labels, name = "Method") +
    theme_minimal(base_size = 14) +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top") +
    ylim(c(-630,0.9)) +
    scale_y_break(c(-400, -30), scale=3) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      axis.title.y.right = element_blank()
    )

  ggsave(result_path_simu("rel_mrae_margin_to_spline_w", NULL, gp_diff_flag, marss_diff_flag, ext = ".jpeg"), plot = p, width = 6, height = 5, dpi = 300)
  
  p <- ggplot(df_stats[df_stats$method %in% subset_methods & df_stats$type=="o",],
                    aes(x = pdt, y = mean_mrae, col = method, shape = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +  # Scatter points
    geom_text(
      data = df_stats[df_stats$method %in% subset_methods & df_stats$type=="o" & df_stats$method=="siap",],
      aes(label = sprintf("%.3f", mean_mrae)),
      vjust = -1, hjust = -0.5,
      size = 3,
      show.legend = FALSE
    ) +
    geom_errorbar(aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae), size = 0.5,
                  width = 0.3, position = position_dodge(width = 0.5)) +  # Error bars
    labs(
      title = "scattered",
      x = "Missingness ratio",
    y = "relative MRAE margin",
      color = "Method"
    ) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_shape_manual(values = method_shapes, labels = method_labels, name = "Method") +
    theme_minimal(base_size = 14) +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top") +
    ylim(c(-30, 1.5))

  ggsave(result_path_simu("rel_mrae_margin_to_spline_o", NULL, gp_diff_flag, marss_diff_flag, ext = ".jpeg"), plot = p, width = 6, height = 5, dpi = 300)
}


type_labels <- c(w = "downtime", o = "scattered")

for (diff_flag in c(T,F)) {
  gp_diff_flag <- diff_flag
  marss_diff_flag <- diff_flag
  df <- get_df_overall_mrae(gp_diff_flag, marss_diff_flag)
  
  df_stats <- df %>%
    group_by(type, method, pdt) %>%
    summarise(
      mean_mrae = mean(rel_margin_mrae, na.rm = T),
      wid_mrae = sd(rel_margin_mrae, na.rm = T) / sqrt(sum(!is.na(rel_margin_mrae))) * qnorm(0.975),
      .groups = "drop"
    )
  
  p_combined <- df_stats[df_stats$method %in% subset_methods, ] %>%
    mutate(type = factor(type, levels = names(type_labels))) %>%
    ggplot(aes(x = pdt, y = mean_mrae, col = method, shape = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(
      aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae),
      width = 0.5, position = position_dodge(width = 0.5), size = 0.5
    ) +
    geom_text(
      data = df_stats[df_stats$method %in% subset_methods & df_stats$type=="w" & df_stats$method=="siap",],
      aes(label = sprintf("%.3f", mean_mrae)),
      vjust = -1, hjust = -.1,
      size = 3,
      show.legend = FALSE
    ) +
    geom_text(
      data = df_stats[df_stats$method %in% subset_methods & df_stats$type=="o" & df_stats$method=="siap",],
      aes(label = sprintf("%.3f", mean_mrae)),
      vjust = -1, hjust = -0.1,
      size = 3,
      show.legend = FALSE
    ) +
    labs(
      # title = "Relative MRAE margin across methods",
      y = "relative MRAE margin",
      x = "Missingness ratio",
      color = "Method"
    ) +
    facet_wrap(~type, nrow = 1, labeller = as_labeller(type_labels)) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_shape_manual(values = method_shapes, labels = method_labels, name = "Method") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top") +
    guides(color = guide_legend(byrow = TRUE, nrow = 1, title.position = "left")) +
    scale_y_break(c(-400, -30), scales = 3) +
    theme(
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      axis.title.y.right = element_blank()
    ) +
    ylim(c(-650, 3))

  ggsave(result_path_simu("rel_mrae_margin_to_spline_combined", NULL, gp_diff_flag, marss_diff_flag, ext = ".jpeg"), plot = p_combined, width = 6, height = 4, dpi = 300)

}

for (diff_flag in c(T,F)) {
  gp_diff_flag <- diff_flag
  marss_diff_flag <- diff_flag
  plot_data_all <- NULL

  for (pdt in c(0.1, 0.3, 0.5)) {
    mrae_stats <- readRDS(file = result_path_simu("mrae_stats", pdt, gp_diff_flag, marss_diff_flag))
    df <- readRDS(file = result_path_simu("overall_rel_mrae_margin", pdt, gp_diff_flag, marss_diff_flag))
    mrae_stats <- mrae_stats %>% 
      left_join(
        df %>% 
        group_by(method) %>%
        summarise(nrepl = n()/2) %>%
          ungroup(), 
        by = c("method")
      )
  
    plot_data <- mrae_stats %>%
      mutate(
        lower_w = mean_mrae.w - 1.96 * sd_mrae.w / sqrt(nrepl),
        upper_w = mean_mrae.w + 1.96 * sd_mrae.w / sqrt(nrepl),
        lower_o = mean_mrae.o - 1.96 * sd_mrae.o / sqrt(nrepl),
        upper_o = mean_mrae.o + 1.96 * sd_mrae.o / sqrt(nrepl)
      )

    plot_data <- plot_data %>%
      mutate(
        pdt = factor(pdt, levels = c(0.1, 0.3, 0.5), labels = c("0.1", "0.3", "0.5"))
      ) %>%
      select(wvl, method, pdt, mean_mrae.w, mean_mrae.o, lower_w, upper_w, lower_o, upper_o) %>%
      pivot_longer(
        cols = c(mean_mrae.w, mean_mrae.o),
        names_to = "type",
        values_to = "mean_mrae"
      ) %>%
      mutate(
        lower = ifelse(type == "mean_mrae.w", lower_w, lower_o),
        upper = ifelse(type == "mean_mrae.w", upper_w, upper_o),
        type = ifelse(type == "mean_mrae.w", "downtime", "scattered")
      )

    plot_data_all <- bind_rows(plot_data_all, plot_data)
  }

  plot_data_all <- plot_data_all %>%
    filter(method %in% subset_methods) %>%
    mutate(type = factor(type, levels = c("downtime", "scattered")))

  p_combined <- ggplot(plot_data_all) +
    #geom_ribbon(aes(x = wvl, ymin = lower, ymax = upper, fill = method), alpha = 0.2) +
    geom_line(aes(x = wvl, y = mean_mrae, color = method, group = method)) +
    facet_grid(
      pdt ~ type
    ) +
    labs(
      x = "Wavelength (nm)",
      y = "Mean Relative Absolute Error",
      color = "Method"
      # fill = "Method"
    ) +
    theme_light(base_size = 14) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_y_log10(sec.axis = dup_axis(name = "Missingness ratio", breaks = NULL, labels = NULL)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE)
      # fill = guide_legend(nrow = 1, byrow = TRUE)
    )

  ggsave(
    result_path_simu("mrae_wvl_combined", NULL, gp_diff_flag, marss_diff_flag, ext = ".jpeg"),
    plot = p_combined, width = 8, height = 6, dpi = 300
  )
}

## ---------- component-wise analysis ------------

type_labels <- c(w = "downtime", o = "scattered")

df <- get_df_overall_mrae(T, T)
  
df_stats <- df %>%
  group_by(type, method, pdt) %>%
  summarise(
    mean_mrae = mean(rel_margin_mrae, na.rm = T),
    wid_mrae = sd(rel_margin_mrae, na.rm = T) / sqrt(sum(!is.na(rel_margin_mrae))) * qnorm(0.975),
    .groups = "drop"
  ) %>% 
  filter(method %in% siap_methods) %>%
  mutate(method = factor(method, levels = c("si_trend_cov", "sia_trend", "s1_si_trend", "si_trend", "si", "sia", "siap")))

p_combined <- df_stats %>%
  mutate(type = factor(type, levels = names(type_labels))) %>%
  ggplot(aes(x = pdt, y = mean_mrae, col = method, shape = method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(
    aes(ymin = mean_mrae - wid_mrae, ymax = mean_mrae + wid_mrae),
    width = 0.5, position = position_dodge(width = 0.5), size = 0.5
  ) +
  geom_text(
    data = df_stats[df_stats$type=="w" & df_stats$method=="siap",],
    aes(label = sprintf("%.3f", mean_mrae)),
    vjust = -1, hjust = -.1,
    size = 3,
    show.legend = FALSE
  ) +
  geom_text(
    data = df_stats[df_stats$type=="o" & df_stats$method=="siap",],
    aes(label = sprintf("%.3f", mean_mrae)),
    vjust = -1, hjust = -0.1,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    y = "relative MRAE margin",
    x = "Missingness ratio",
    color = "Method"
  ) +
  facet_wrap(~type, nrow = 1, labeller = as_labeller(type_labels)) +
  scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
  scale_shape_manual(values = method_shapes, labels = method_labels, name = "Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top") +
  guides(color = guide_legend(byrow = FALSE, nrow = 2, title.position = "left")) +
  theme(
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.title.y.right = element_blank()
  ) + 
  ylim(c(-3.5,1.5))

ggsave(result_path_simu("rel_mrae_margin_to_spline_siap", ext = ".jpeg"), plot = p_combined, width = 7, height = 4, dpi = 300)


## ---------- UQ --------------
ex <- apply(xt, 1, mean)
wvl <- as.numeric(rownames(xt))
cap_val <- 0.1 # Cap MARSS spikes only for visualization (keep raw values in `rel_length`).

for (diff_flag in c(T,F)) {
  gp_diff_flag <- diff_flag
  marss_diff_flag <- diff_flag
  
  combined_data <- NULL
  for (pdt in c(0.1, 0.3, 0.5)) {
    plot_data <- readRDS(file = result_path_simu("length_stats", pdt, gp_diff_flag, marss_diff_flag)) %>%
      filter(method %in% subset_methods)
    
    plot_data_w <- plot_data %>%
      mutate(
        pdt = factor(pdt, levels = c(0.1, 0.3, 0.5), labels = c("0.1", "0.3", "0.5")),
        missingness = "Downtime",
        rel_length = mean_length.w / ex
      )
    
    plot_data_o <- plot_data %>%
      mutate(
        pdt = factor(pdt, levels = c(0.1, 0.3, 0.5), labels = c("0.1", "0.3", "0.5")),
        missingness = "Scattered",
        rel_length = mean_length.o / ex
      )
    
    combined_data <- bind_rows(combined_data, plot_data_w, plot_data_o)
  }
  
  
  combined_data <- combined_data %>%
    mutate(
      capped = method == "marss" & rel_length > cap_val,
      rel_length_plot = ifelse(capped, cap_val, rel_length)
    )
  
  p_combined <- ggplot(combined_data) +
    geom_line(aes(x = wvl, y = rel_length_plot, color = method, group = method)) +
    geom_point(
      data = combined_data %>% filter(capped),
      aes(x = wvl, y = rel_length_plot, color = "gray"),
      shape = 4, size = 1.8, stroke = 0.9
    ) +
    facet_grid(
      pdt ~ missingness
    ) +
    labs(
      x = "Wavelength (nm)",
      y = "Relative Interval Width",
      color = "Method"
    ) +
    theme_light(base_size = 14) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_y_log10(
      sec.axis = dup_axis(name = "Missingness ratio", breaks = NULL, labels = NULL)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
  
  ggsave(
    result_path_simu("rel_length_wvl_combined", NULL, gp_diff_flag, marss_diff_flag, ext = ".jpeg"),
    plot = p_combined, width = 8, height = 6, dpi = 300
  )
}

for (diff_flag in c(T, F)) {
  gp_diff_flag <- diff_flag
  marss_diff_flag <- diff_flag
  df_combined <- NULL

  for (pdt in c(0.1, 0.3, 0.5)) {
    df <- readRDS(file = result_path_simu("coverage_wvl_stats", pdt, gp_diff_flag, marss_diff_flag)) %>%
      filter(method %in% subset_methods)
    df$wvl <- as.numeric(df$wvl)

    df0 <- readRDS(file = result_path_simu("overall_rel_mrae_margin", pdt, gp_diff_flag, marss_diff_flag))
    df <- df %>%
      left_join(
        df0 %>%
          group_by(method) %>%
          summarise(nrepl = n() / 2, .groups = "drop"),
        by = c("method")
      ) %>%
      mutate(
        lower_w = w - 1.96 * sd_w / sqrt(nrepl),
        upper_w = w + 1.96 * sd_w / sqrt(nrepl),
        lower_o = o - 1.96 * sd_o / sqrt(nrepl),
        upper_o = o + 1.96 * sd_o / sqrt(nrepl)
      )

    df_combined <- rbind(df_combined, mutate(df, pdt = pdt))
  }

  plot_df <- df_combined %>%
    select(wvl, method, pdt, w, o, lower_w, upper_w, lower_o, upper_o) %>%
    pivot_longer(cols = c(w, o), names_to = "missing_type", values_to = "coverage") %>%
    mutate(
      lower = if_else(missing_type == "w", lower_w, lower_o),
      upper = if_else(missing_type == "w", upper_w, upper_o),
      missing_type = factor(
        missing_type,
        levels = c("w", "o"),
        labels = c("downtime", "scattered")
      ),
      pdt = factor(
        pdt,
        levels = c(0.1, 0.3, 0.5),
        labels = c("0.1", "0.3", "0.5")
      ),
      method_group = if_else(method == "marss", "MARSS", "GP, SI, SIAP"),
      method_group = factor(method_group, levels = c("GP, SI, SIAP", "MARSS"))
    )

  combined_plot <- ggplot(plot_df, aes(x = wvl, y = coverage, color = method, group = method)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.2, color = NA) +
    geom_line() +
    geom_hline(
      data = plot_df %>%
        dplyr::filter(method_group == "GP, SI, SIAP") %>%
        dplyr::distinct(pdt, method_group, missing_type) %>%
        dplyr::mutate(yint = 0.95),
      aes(yintercept = yint),
      linetype = "dashed",
      color = "darkgrey",
      linewidth = 1.2,
      alpha = 0.8
    ) +
    labs(x = "Wavelength (nm)", y = expression(frac(1, N[repl]) * sum(AveCov["i" * "," * "\u2113"], "\u2113"==1, N[repl]))) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(sec.axis = dup_axis(name = "Missingness ratio", breaks = NULL, labels = NULL)) +
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_fill_manual(values = method_colors, labels = method_labels, name = "Method") +
    theme_light(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal"
    ) +
    ggh4x::facet_nested(
      rows = vars(pdt, method_group),
      cols = vars(missing_type),
      scales = "free_y",
      labeller = labeller(method_group = c("GP, SI, SIAP" = "", "MARSS" = "")),
      nest_line = element_blank(),
      strip = ggh4x::strip_nested(
        size = "variable",
        by_layer_y = TRUE,
        text_y = list(element_text(), element_blank()),
        background_y = list(element_rect(), element_blank())
      )
    )

  ggsave(
    result_path_simu("coverage_wvl", NULL, gp_diff_flag, marss_diff_flag, ext = ".jpeg"),
    plot = combined_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}



# -------------- Table -------------
for (diff_flag in c(T, F)) {
  gp_diff_flag <- diff_flag
  marss_diff_flag <- diff_flag
  coverage_tb <-NULL
  
  for (pdt in c(0.1,0.3,0.5)) {
    coverage_stats <- readRDS(file = result_path_simu("coverage_stats", pdt, gp_diff_flag, marss_diff_flag))
    df0 <- readRDS(file = result_path_simu("overall_rel_mrae_margin", pdt, gp_diff_flag, marss_diff_flag))
    coverage_stats <- coverage_stats %>%
      left_join(
        df0 %>%
          group_by(method) %>%
          summarise(nrepl = n() / 2, .groups = "drop"),
        by = c("method")
      ) %>% 
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
  # print(coverage_summary)
  write.table(coverage_summary, file = result_path_simu("coverage_table", NULL, gp_diff_flag, marss_diff_flag, ext = ".csv"), row.names = F)
}
