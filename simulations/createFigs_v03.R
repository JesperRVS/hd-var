# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggplot2", "reshape")
lapply(libplt, require, character.only = TRUE)
# devtools::install_github("stefano-meschiari/latex2exp")
library("latex2exp")

# Load workspace
# load("simulations_workspace_10_MC_100_to_500_n_4_to_16_p_3_h.RData")
load("simulations_workspace_80_MC_100_to_500_n_16_to_128_p_3_h_LINUX_with_num_upd.RData")

## == HELPER FUNCTIONS ##
# Function providing locations of methods to plot
which_methods <- function(met_plt) {
  nummet_plt <- length(met_plt)
  where_met_plt <- numeric(nummet_plt)
  for (thismet in 1:nummet_plt) {
    where_met_plt[thismet] <- which(methods == met_plt[thismet])
  }
  return(where_met_plt)
}

# Create labels (via functions) for plots
plab <- function(string) {
  sprintf("p = %s", string) # note: labeller expects a string - not integer
}

# Specify statistics to plot and methods to include
stats_plt <- c("max", "mean", "median", "sd") # statistics to plot (fed to apply)

# Function providing string for statistic (included in y-label)
stat_string <- function(stat) {
  switch(stat,
         "max" = "Maximum",
         "mean" = "Average",
         "median" = "Median",
         "var" = "Variance",
         "sd" = "Standard Deviation",
         stop("Invalid statistic"))
}

# Function taking statistic and methods to plot, returning ggplot
t_str <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"
plot_stat_p_by_design <- function(stat, met_plt) {
  where_met_plt <- which_methods(met_plt)
  stat_max_ell2_errors <- apply(max_ell2_errors, c(2, 3, 4, 5), stat)
  df <- reshape::melt(stat_max_ell2_errors[, , , method = where_met_plt])
  plt <- ggplot(df, aes(x = n, y = value, color = as.factor(method))) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(design ~ p, labeller = labeller(p = plab),
               scales = "free_y") +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str))
    ) +
    # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
            legend.key = element_rect(colour = "transparent", fill = NA),
            legend.box.background = element_rect(colour = "black"),
            legend.title = element_blank(),
            legend.text = element_text(size = rel(0.5)),
            legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
            legend.spacing = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            panel.spacing.y = unit(3, "mm"),
            plot.title = element_text(size = rel(0.5)),
            plot.caption = element_text(hjust = 0, size = rel(0.5))
    )
  return(plt)
}

# 1. Comparisons w/o refitting (i.e. no "Post" in name)
met_plt <- c("Lasso", "SqrtLasso")
# met_plt <- c("Lasso", "AICLasso", "BICLasso", "SqrtLasso")
# met_plt <- methods[!grepl("Post", methods)]

for (thisstat in seq_along(stats_plt)) {
  plt <- plot_stat_p_by_design(stats_plt[thisstat], met_plt)
  file_name <- paste("plain_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}

# 2. Comparisons w/ refitting (i.e. "Post" in name)
met_plt <- c("PostLasso", "PostBICLasso", "PostSqrtLasso")
# met_plt <- c("PostLasso", "PostAICLasso", "PostBICLasso", "PostSqrtLasso")
# met_plt <- methods[grepl("Post", methods)]
for (thisstat in seq_along(stats_plt)) {
  plt <- plot_stat_p_by_design(stats_plt[thisstat], met_plt)
  file_name <- paste("post_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_p_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}

# Function taking statistic and methods to plot, returning ggplot
# with design as the row facet, method as the column facet, n as x-axis
# value as y-axis, and color by p
t_str <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"
plot_stat_method_by_design <- function(stat, met_plt) {
  where_met_plt <- which_methods(met_plt)
  stat_max_ell2_errors <- apply(max_ell2_errors, c(2, 3, 4, 5), stat)
  df <- reshape::melt(stat_max_ell2_errors[, , , method = where_met_plt])
  g <- guide_legend(title = TeX("$p =$"))
  plt <- ggplot(df, aes(x = n, y = value, color = as.factor(p))) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    facet_grid(design ~ method, labeller = labeller(p = plab),
               scales = "free_y", breaks = scales::pretty_breaks(n = 5)) +
    labs(
      x = TeX(r"($n$)"),
      y = TeX(paste0("Monte Carlo ", stat_string(stat), t_str))
    ) +
    # set number of yticks to 5
    # scale_y_continuous() +
    theme_bw() +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),
            legend.key = element_rect(colour = "transparent", fill = NA),
            legend.box.background = element_rect(colour = "black"),
            legend.title = element_text(size = rel(0.5), hjust = 0.5),
            legend.text = element_text(size = rel(0.5)),
            legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
            legend.spacing = unit(0, "mm"),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"),
            panel.spacing.y = unit(3, "mm"),
            plot.title = element_text(size = rel(0.5)),
            plot.caption = element_text(hjust = 0, size = rel(0.5))
    ) +
    guides(color = g)
  return(plt)
}

# 1. Comparisons w/o refitting (i.e. no "Post" in name)
met_plt <- c("Lasso", "SqrtLasso")

for (thisstat in seq_along(stats_plt)) {
  plt <- plot_stat_method_by_design(stats_plt[thisstat], met_plt)
  file_name <- paste("plain_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_method_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}

# 2. Comparisons all Lasso and SqrtLasso methods
met_plt <- c("Lasso", "SqrtLasso", "PostLasso", "PostSqrtLasso")

for (thisstat in seq_along(stats_plt)) {
  plt <- plot_stat_method_by_design(stats_plt[thisstat], met_plt)
  file_name <- paste("all_lasso_sqrtlasso_methods_",
                     stats_plt[thisstat], "_max_ell2_errors_method_by_design_",
                     nummc, "_MC_", min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 8, height = 12)
}


# OLD BELOW THIS LINE
# df_mean <- reshape::melt(mean_max_ell2_errors)
# df_sd <- reshape::melt(sd_max_ell2_errors)
# names(df_mean)[names(df_mean) == "value"] <- "mean"
# df <- cbind(df_mean, df_sd["value"])
# names(df)[names(df) == "value"] <- "sd"

# # Function taking methods to plot, returning ggplot of means
# # including error bars for plus/minus one standard deviation
# t_str <- " of $\\max_{i\\in[p]}\\|\\widehat{\\beta}_i-\\beta_{{0}{i}}\\|_{\\ell_2}$"
# plot_mean_p_by_design_with_error_bars <- function(met_plt) {
#   where_met_plt <- which_methods(met_plt)
#   these_max_ell2_errors <- max_ell2_errors[, , , , method = where_met_plt]
#   mean_max_ell2_errors <- apply(these_max_ell2_errors, c(2, 3, 4, 5), mean)
#   sd_max_ell2_errors <- apply(these_max_ell2_errors, c(2, 3, 4, 5), sd)
#   df_mean <- reshape::melt(mean_max_ell2_errors)
#   df_sd <- reshape::melt(sd_max_ell2_errors)
#   names(df_mean)[names(df_mean) == "value"] <- "mean"
#   df <- cbind(df_mean, df_sd["value"])
#   names(df)[names(df) == "value"] <- "sd"
#   plt <- ggplot(df, aes(x = n, y = mean, color = as.factor(method))) +
#     geom_line(linewidth = 0.5) +
#     geom_point(size = 2) +
#     geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
#     facet_grid(design ~ p, labeller = labeller(p = plab)) +
#     labs(
#       x = TeX(r"($n$)"),
#       y = TeX(paste0("Monte Carlo mean and SD", t_str))
#     ) +
#     theme_bw() +
#     theme(legend.justification = c(1, 1), legend.position = c(1, 1),
#           legend.key = element_rect(colour = "transparent", fill = NA),
#           legend.box.background = element_rect(colour = "black"),
#           legend.title = element_blank(),
#           legend.text = element_text(size = rel(0.5)),
#           legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
#           legend.spacing = unit(0, "mm"),
#           legend.spacing.x = unit(0, "mm"),
#           legend.spacing.y = unit(0, "mm"),
#           panel.spacing.y = unit(3, "mm"),
#           plot.title = element_text(size = rel(0.5)),
#           plot.caption = element_text(hjust = 0, size = rel(0.5))
#     )
# }

# # 1. Comparisons w/o refitting (i.e. no "Post" in name)
# met_plt <- c("Lasso", "BICLasso", "SqrtLasso")
# plt_plain <- plot_mean_p_by_design_with_error_bars(met_plt)
# file_name <- paste("plain_methods_mean_with_error_bars_max_ell2_errors_p_by_design_",
#                    nummc, "_MC_", min(nvec), "_to_", max(nvec),
#                    "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
# ggsave(plt_plain, filename = paste0("img/", file_name, ".png"),
#        width = 8, height = 12)

# # Function like the one above but without error bars
# plot_mean_p_by_design <- function(met_plt) {
#   where_met_plt <- which_methods(met_plt)
#   these_max_ell2_errors <- max_ell2_errors[, , , , method = where_met_plt]
#   mean_max_ell2_errors <- apply(these_max_ell2_errors, c(2, 3, 4, 5), mean)
#   df <- reshape::melt(mean_max_ell2_errors)
#   plt <- ggplot(df, aes(x = n, y = value, color = as.factor(method))) +
#     geom_line(linewidth = 0.5) +
#     geom_point(size = 2) +
#     facet_grid(design ~ p, labeller = labeller(p = plab)) +
#     labs(
#       x = TeX(r"($n$)"),
#       y = TeX(paste0("Monte Carlo mean", t_str))
#     ) +
#     theme_bw() +
#     theme(legend.justification = c(1, 1), legend.position = c(1, 1),
#           legend.key = element_rect(colour = "transparent", fill = NA),
#           legend.box.background = element_rect(colour = "black"),
#           legend.title = element_blank(),
#           legend.text = element_text(size = rel(0.5)),
#           legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
#           legend.spacing = unit(0, "mm"),
#           legend.spacing.x = unit(0, "mm"),
#           legend.spacing.y = unit(0, "mm"),
#           panel.spacing.y = unit(3, "mm"),
#           plot.title = element_text(size = rel(0.5)),
#           plot.caption = element_text(hjust = 0, size = rel(0.5))
#     )
# }

# # 1. Comparisons w/o refitting (i.e. no "Post" in name)
# met_plt <- c("Lasso", "BICLasso", "SqrtLasso")
# plt_plain <- plot_mean_p_by_design(met_plt)
# file_name <- paste("plain_methods_mean_max_ell2_errors_p_by_design_",
#                    nummc, "_MC_", min(nvec), "_to_", max(nvec),
#                    "_n_", min(pvec), "_to_", max(pvec), "_p", sep = "")
# ggsave(plt_plain, filename = paste0("img/", file_name, ".png"),
#        width = 8, height = 12)