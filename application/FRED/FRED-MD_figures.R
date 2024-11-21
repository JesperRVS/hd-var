# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)

# Load data
load("application_workspace_N_120_qmax_12.RData")

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

met_lab <- function(string) {
  met_list <- c("Lasso" = "Weighted Lasso",
                "AICLasso" = "AIC-Lasso",
                "BICLasso" = "BIC-Lasso",
                "HQICLasso" = "HQIC-Lasso",
                "SqrtLasso" = "Sqrt-Lasso",
                "PostLasso" = "Post-Lasso",
                "PostAICLasso" = "Post-AIC-Lasso",
                "PostBICLasso" = "Post-BIC-Lasso",
                "PostHQICLasso" = "Post-HQIC-Lasso",
                "PostSqrtLasso" = "Post-Sqrt-Lasso")
  return(met_list[string])
}

# Function providing string for statistic (included in y-label)
stat_string <- function(stat) {
  switch(stat,
         "max" = "Maximum",
         "mean" = "Average",
         "median" = "Median",
         "q90" = "90th Quantile",
         "q95" = "95th Quantile",
         "var" = "Variance",
         "sd" = "Standard Deviation",
         stop("Invalid statistic"))
}

# Function providing (function for) statistic
stat_fctn <- function(stat) {
  switch(stat,
         "max" = max,
         "mean" = mean,
         "median" = median,
         "q90" = function(x) quantile(x, probs = 0.9),
         "q95" = function(x) quantile(x, probs = 0.95),
         "var" = var,
         "sd" = sd,
         stop("Invalid statistic"))
}

# Color-blind friendly color palette
cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Shape order
shape_order <- c(0, 1, 2, 5, 4, 3)
# 0=square, 1=circle, 2=triangle, 5=diamond, 4=cross, 3=plus

# Linetype order
linetype_order <- c("solid", "dashed", "dotdash",
                    "dotted", "longdash", "twodash")

# Plotting function
plot_rel_stat_forecast_error <- function(stat, met_plt) {
  where_met_plt <- which_methods(met_plt)
  dimnames(ivwsfes)[["method"]] <-
    met_lab(dimnames(ivwsfes)[["method"]])
  stat_ivwsfes <- apply(ivwsfes, c(2, 3), stat_fctn(stat)) # mean IVWSFE
  rel_stat_ivwsfes <- 100 * stat_ivwsfes / stat_ivwsfes[1, 1] # / VAR(1) LASSO
  rel_stat_iwvsfes_plt <- rel_stat_ivwsfes[, where_met_plt] # methods to plot
  df <- reshape::melt(rel_stat_iwvsfes_plt)
  g <- guide_legend(nrow = 1)
  plt <- ggplot(df, aes(x = q, y = value)) +
    geom_hline(yintercept = 100, linetype = "dotted", color = "black") +
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(shape = method, color = method), size = 2) +
    scale_shape(solid = FALSE) +
    labs(
      x = TeX("Lag length ($q$)"),
      y = TeX(paste(stat_string(stat),
                    "IVWSFE in % of VAR$(1)$ Lasso")),
      ) +
      theme_bw() +
      theme(#legend.justification = c(.99, .01),
        #legend.position = c(.99, .01),
        legend.position = "bottom",
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(.75)),
        legend.margin = margin(c(0, 0, 0, 0), unit = "mm"),
        legend.spacing = unit(0, "mm"),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        panel.spacing.y = unit(3, "mm"),
        plot.title = element_text(size = rel(0.5)),
        plot.caption = element_text(hjust = 0, size = rel(0.5))
      ) +
      scale_x_continuous(breaks = seq(2, numlag, 2), limits = c(1, 12)) +
      scale_colour_manual(values = cb_palette) +
      scale_shape_manual(values = shape_order) +
      scale_linetype_manual(values = linetype_order) +
      guides(color = g, shape = g, linetype = g)
    return(plt)
}

# Plot the mean IVWSFE
met_plt <- c("Lasso", "PostLasso", "SqrtLasso")#, "BICLasso")
stats_plt <- c("mean", "median", "q90", "q95") # statistics to plot
for (thisstat in seq_along(stats_plt)) {
  which_stat <- stats_plt[thisstat]
  plt <- plot_rel_stat_forecast_error(stats_plt[thisstat], met_plt)
  file_name <- paste0("rel_", stats_plt[thisstat],
                      "_ivwsfe_N_", numfore, "_qmax_", numlag)
  ggsave(plt, filename = paste0("img/", file_name, ".png"),
         width = 4, height = 4)
}