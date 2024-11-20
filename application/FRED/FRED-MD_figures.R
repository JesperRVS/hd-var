# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for plotting
libplt <- c("devtools", "ggh4x", "ggplot2", "latex2exp", "reshape", "scales")
lapply(libplt, require, character.only = TRUE)

# Load data
load("application_workspace_N_10_qmax_4.RData")

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
plot_rel_mean_forecast_error <- function(met_plt) {
  where_met_plt <- which_methods(met_plt)
  dimnames(ivwsfes)[["method"]] <-
    met_lab(dimnames(ivwsfes)[["method"]])
  mivwsfes <- apply(ivwsfes, c(2, 3), mean) # mean IVWSFE
  rel_mivwsfes <- mivwsfes / mivwsfes[1, 1] # rel. to VAR(1) LASSO
  rel_miwvsfes_plt <- rel_mivwsfes[, where_met_plt] # methods to plot
  df <- reshape::melt(rel_miwvsfes_plt)
  g <- guide_legend(nrow = 1)
  plt <- ggplot(df, aes(x = q, y = value)) + 
    geom_line(aes(linetype = method, color = method), linewidth = 0.5) +
    geom_point(aes(shape = method, color = method), size = 2) +
    scale_shape(solid = FALSE) +
    labs(
      x = TeX("Lag length ($q$)"),
      y = "MIVWSFE relative to VAR(1) Lasso (= 1)",
    ) +
    theme_bw() +
    theme(
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
    scale_colour_manual(values = cb_palette) +
    scale_shape_manual(values = shape_order) +
    scale_linetype_manual(values = linetype_order) +
    guides(color = g, shape = g, linetype = g)
  return(plt)
}

# Plot the mean IVWSFE
met_plt <- c("Lasso", "PostLasso", "SqrtLasso", "BICLasso")
plt <- plot_rel_mean_forecast_error(met_plt)
file_name <- paste0("rel_mean_ivwsfe_N_", numfore, "_qmax_", numlag)
ggsave(plt, filename = paste0("img/", file_name, ".png"), width = 6, height = 4)